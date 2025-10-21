import os 
import json 
import hashlib 
from pathlib import Path
from sseqs import msa

from fastapi import Request, Form, FastAPI
from fastapi.responses import Response, StreamingResponse, FileResponse
from fastapi.middleware.cors import CORSMiddleware
import uvicorn
import io, tarfile

from fastapi.responses import StreamingResponse
import json
from pydantic import BaseModel
from typing import List, Literal

import argparse
parser = argparse.ArgumentParser(description="Run MSA Server")
parser.add_argument( "-port", '-p', default=8000, type=int, help="Port for server. ")
args, _ = parser.parse_known_args()

app = FastAPI()
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_credentials=True, allow_methods=["*"], allow_headers=["*"],)

chunks = int(os.environ.get('CHUNKS', 4))

# code below supports `boltz predict --msa_server_url=0.0.0.0:8000 --use_msa_server`
result_store: dict[str, Path] = {}
@app.post("/ticket/msa", include_in_schema=False)
@app.post("/ticket/pair", include_in_schema=False)
async def ticket_msa(request: Request, q: str = Form(...)):
    q = str(q).split('\n')

    # grab parameters from environment variables. 
    msas_num = int(os.getenv("MSAS_NUM", 4096))
    msas_fft_rank = int(os.getenv("MSAS_FFT_RANK", 4))
    msas_top_fft = int(os.getenv("MSAS_TOP_FFT", 200))
    msas_top_sw = int(os.getenv("MSAS_TOP_SW", 10))
    msas_top_sw_affine = int(os.getenv("MSAS_TOP_SW_AFFINE", 2))

    # cache results in folder specific to options. 
    cache_folder = f"{chunks}_{msas_num}_{msas_top_fft}_{msas_fft_rank}_{msas_top_sw}_{msas_top_sw_affine}"
    os.makedirs(f"cache_msa/{cache_folder}/", exist_ok=True)

    # fetch all proteins from input 
    heads, proteins = q[::2], q[1::2]
    out_paths = []
    for protein in proteins: 
        protein_hash = hashlib.sha256(protein.encode('utf-8')).hexdigest()
        out_paths.append(f"cache_msa/{cache_folder}/{protein_hash}.a3m")

    # compute unique proteins and output paths 
    unique_proteins, a3m_paths = [], []
    for protein in list(set(proteins)):
        protein_hash = hashlib.sha256(protein.encode('utf-8')).hexdigest()
        a3m_path = f"cache_msa/{cache_folder}/{protein_hash}.a3m" 
        if os.path.exists(a3m_path): continue
        unique_proteins.append(protein)
        a3m_paths.append(a3m_path)

    # compute MSA for all proteins 
    if len(unique_proteins)>0: 
        msa(unique_proteins, 
                a3m_paths,
                fft_rank=msas_fft_rank,
                top_fft=msas_top_fft, 
                top_sw=msas_top_sw, 
                top_sw_affine=msas_top_sw_affine, 
                num_msas=msas_num,
                bs=100_000_000, )

    # fix header in output format. 
    for i, a3m_path in enumerate(out_paths): 
        a = open(a3m_path, 'r').read().replace("@",'X')
        lines = a.split('\n')
        lines[0] = q[i*2] 
        open(a3m_path, 'w').write("\n".join(lines))

    # add seperator (?null ascii character?). 
    uniref_content = b"\x00\n".join(Path(p).read_bytes() for p in out_paths)

    # satisfy the output format
    tar_id = protein_hash 
    tar_path = Path(tar_id.replace('.a3m','_tar'))
    tar_path.mkdir(exist_ok=True)
    tar_file = tar_path / f"{tar_id}.tar.gz"

    with tarfile.open(tar_file, "w:gz") as tar:
        is_pair = request.url.path.endswith("/ticket/pair")
        name = "pair.a3m" if is_pair else "uniref.a3m"
        ti = tarfile.TarInfo(name=name)
        ti.size = len(uniref_content)
        tar.addfile(ti, io.BytesIO(uniref_content))
        env_name = "bfd.mgnify30.metaeuk30.smag30.a3m"
        ti2 = tarfile.TarInfo(name=env_name)
        ti2.size = 0
        tar.addfile(ti2, io.BytesIO(b""))

    result_store[tar_id] = tar_file

    return {"status": "COMPLETE", "id": tar_id}

# Endpoint to serve the tar.gz back to run_mmseqs2
@app.get("/result/download/{tar_id}", include_in_schema=False)
async def result_download(tar_id: str):
    path = result_store.get(tar_id)
    if path is None or not path.exists():
        return Response(status_code=404)
    return FileResponse(path, media_type="application/gzip", filename="out.tar.gz")



class MSARequest(BaseModel): 
    from pydantic import Field
    sequences: List[str] = Field(["ADAM", "ADA"], \
                                 description="List of input sequences, will compute one Multiple Sequence Alignment for each input.")
    msas_num: int = Field(4096,     description="Number of lines in output MSA, at most 4096. ")
    fft_rank: int = Field(1,        description="Iterations of initial FFT filter (informally, rank of SVD(blosum62*query_sequence)). ")
    top_fft: int = Field(10,        description="Stage 1:  Pass on `1//top_fft` sequences from FFT filter to gapped Smith-Waterman gap_cost=-1. ")
    top_sw: int = Field(2,          description="Stage 2:  Pass on `1//top_sw` sequences from Smith-Waterman gap_cost=-1 to gapped Smith-Waterman with open_gap=-11 and extend_gap=-1. ")
    top_sw_affine: int = Field(2,   description="Stage 3:  Pass on `1//top_sw_affine` from Smith-Waterman with open_gap=-11 and extend_gap=-1 sequences to build .a3m file. ")

    dataset: Literal["uniref30_2302+bfd_mgy_colabfold"] = \
    Field(..., description="Dataset to use for alignment, request new datasets by email: Alexander.Mathiasen@gmail.com")

class MSAResult(BaseModel): 
    alignment: List[str]
class ProgressEvent(BaseModel):
    progress: float
    message: str

@app.post("/msa", response_model=MSAResult)
def get_msa(req: MSARequest):
    """ Returns Multiple Sequence Alignment for the given input sequences.  """

    cache_folder = f"{chunks}_{req.msas_num}_{req.top_fft}_{req.fft_rank}_{req.top_sw}_{req.top_sw_affine}"
    os.makedirs(cache_folder, exist_ok=True)
    proteins = req.sequences
    filenames_a3m = []
    for protein in proteins:
        protein_hash = hashlib.sha256(protein.encode('utf-8')).hexdigest()
        filenames_a3m.append(f"{cache_folder}/{protein_hash}.a3m")

    def event_stream():
        try: 
            msa(req.sequences, filenames_a3m=filenames_a3m, 
                            num_msas=req.msas_num, 
                            bs=100_000_000,
                            fft_rank=req.fft_rank, 
                            top_fft=req.top_fft,  
                            top_sw=req.top_sw,  
                            top_sw_affine=req.top_sw_affine)

            dct = {}
            for i, file in enumerate(filenames_a3m):
                msa_string = open(file, "r").read()
                dct[i] = msa_string

            yield json.dumps({
                'type': 'data_msa', 
                'data': dct, 
                'num_msas': len(proteins),
                'cache_folder': cache_folder
            }) 
        except Exception as e:
            yield json.dumps({"error": f"MSA did not succeed: {e}."})

    return StreamingResponse(event_stream(), media_type="application/x-ndjson")

uvicorn.run(app, host="0.0.0.0", port=args.port, reload=False, access_log=False, workers=1)