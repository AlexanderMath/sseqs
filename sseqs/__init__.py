# Lazy imports to avoid CUDA compilation on every import
def __getattr__(name):
    if name == 'sw':
        from .sw import sw
        return sw
    if name == 'sw_profile_cuda':
        from .sw import sw_profile_cuda
        return sw_profile_cuda
    if name == 'pssm_from_file':
        from .utils import pssm_from_file
        return pssm_from_file
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
