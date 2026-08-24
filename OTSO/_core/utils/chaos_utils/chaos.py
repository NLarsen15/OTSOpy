from chaos_method import get_chaos_g_h
from chaos_download import download_mat_only
from datetime import datetime

CHAOSDOI = "10.5281/zenodo.13950012"

dt = datetime(2013, 9, 1, 0, 0, 0)
files = download_mat_only(CHAOSDOI)
g, h = get_chaos_g_h(dt)

print(g)
print(h)