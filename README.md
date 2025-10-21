<img src='images/goal.png'/>
MSA by ML researchers for ML researchers --️ in pytorch/cuda ❤️

```
# install
pip install sseqs
wget https://foldify.org/uniref_bfd_mgy_cf.xbit 

# python 
from sseqs import msa
msa("HPETLVKVKDAEDQLGARVG"*10, "msa.a3m")
db_len=998 q_len=200: 100%|█████████| 121/121 [00:11<00:00, 10.40GB/s]

# boltz2 msa-server 🔥
DBPATH=uniref_bfd_mgy_cf.xbit python server.py --port 8000
boltz predict demo.fasta --msa_server_url http://0.0.0.0:8000 --use_msa_server
```

No need for <a href="https://instances.vantage.sh/aws/ec2/x2gd.16xlarge?currency=USD">$5h/h</a> server with 1000GB RAM. 
Developed for <a href="https://cloud.vast.ai/">$0.3/h</a> rtx4090+128GB RAM. 

# limitations
- no protein pairing 
- sequence length < 1000 (working on 2048)
- 128gb RAM (working on 64GB w/ compression) 
- no <16GB approximate version (yet)
- no evals (working on runs'n'poses + antibody)
