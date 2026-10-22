###############################################################################
# Export real trajectories and their CVAE reconstructions for the appendix
# figure of the paper (reconstruction error, 48 trajectories).
#
# The preprocessing is the one of train_traj_vae.py (weekdays only, projection
# in meters, scale = 1e4). The reconstruction is the deterministic output of
# the decoder evaluated at the posterior mean: G(mu_m, c_m).
#
# Usage (from the py_workflow folder):
#   python export_reconstruction.py --nz 3 --n 48
# Output: ../R/pydata/reconstruction_cvae<nz>.csv
###############################################################################

import argparse
from datetime import datetime

import numpy as np
import pandas as pd
import torch

import cvae_module

parser = argparse.ArgumentParser(description="Export CVAE reconstructions")
parser.add_argument("--data_path", default="../R/Rdata/xgv_100_14-18.csv")
parser.add_argument("--save_folder", default="save")
parser.add_argument("--out_folder", default="../R/pydata")
parser.add_argument("--nz", default=3, type=int)
parser.add_argument("--n", default=48, type=int)
parser.add_argument("--seed", default=1234, type=int)
args = parser.parse_args()

###############################################################################
# Treating the data (same as train_traj_vae.py)
###############################################################################
data = pd.read_csv(args.data_path)
padding = 24
marees = data.maree.unique()

traj_all = np.zeros((len(marees), 2, padding))
day_id = np.zeros(len(marees))
day = np.zeros(len(marees))
for i, (tt, traj) in enumerate(data.groupby("maree", sort=False)):
    xy = traj.loc[:, ("x", "y")].to_numpy()
    traj_all[i, :, 0:xy.shape[0]] = xy.T
    day_id[i] = datetime.weekday(datetime.fromisoformat(traj.time.iloc[0]))
    day[i] = datetime.strptime(traj.time.iloc[1], "%Y-%m-%d").timetuple().tm_yday

keep = day_id < 5
marees = marees[keep]
day = day[keep]

R = 6377726  # earth radius
scale = 10000  # empirical choice, as in training
traj_std = np.zeros((keep.sum(), 2, padding))
traj_std[:, 0, :] = np.pi * R * (traj_all[keep, 0, :] + 4.28) / 180 / scale
traj_std[:, 1, :] = np.pi * R * (traj_all[keep, 1, :] - 47.78) / 180 / scale

c_day = torch.tensor(
    np.stack([np.cos(day / 365.25 * 2 * np.pi), np.sin(day / 365.25 * 2 * np.pi)], 1),
    dtype=torch.float)[:, :, None]
real = torch.tensor(traj_std, dtype=torch.float)

###############################################################################
# Loading the CVAE and reconstructing the trajectories from the posterior mean
###############################################################################
device = "cpu"
E = cvae_module.cdepthCNNEncoder(nz=args.nz, nc=2).to(device)
G = cvae_module.cCNNGenerator(nz=args.nz, nc=2).to(device)
E.load_state_dict(torch.load("%s/E%d_cvae_std_10k_adam_epoch_1000.pth" % (args.save_folder, args.nz),
                             map_location=torch.device(device)))
G.load_state_dict(torch.load("%s/G%d_cvae_std_10k_adam_epoch_1000.pth" % (args.save_folder, args.nz),
                             map_location=torch.device(device)))
E.eval()
G.eval()
with torch.no_grad():
    _, mu, _ = E(real, c_day)
    x_hat = G(mu[:, :, None], c_day)

###############################################################################
# Export a random sample of trajectories (long format)
###############################################################################
rng = np.random.default_rng(args.seed)
ind = rng.choice(real.shape[0], size=args.n, replace=False)

def to_long(arr, kind):
    return pd.DataFrame({
        "id": np.repeat(np.arange(1, args.n + 1), padding),
        "maree": np.repeat(marees[ind], padding),
        "h": np.tile(np.arange(1, padding + 1), args.n),
        "x": arr[ind, 0, :].reshape(-1),
        "y": arr[ind, 1, :].reshape(-1),
        "type": kind})

out = pd.concat([to_long(real.numpy(), "real"), to_long(x_hat.numpy(), "reconstructed")])
out_file = "%s/reconstruction_cvae%d.csv" % (args.out_folder, args.nz)
out.to_csv(out_file, index=False)
print("Saved", out_file)
