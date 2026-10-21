import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import torch
from torch.utils.data import DataLoader, Dataset
from torch import nn
from torch import optim
import torch.nn.functional as F
from datetime import datetime
import cond_vae



np.random.seed(1234)
torch.manual_seed(1234)
torch.cuda.manual_seed(1234)


###############################################################################
# Parsing script arguments
###############################################################################
data_file = "~/Documents/maccoloc/maccoRoc/data/xgv/xgv_100_14-18.csv"

parser = argparse.ArgumentParser(description="cVAE")
parser.add_argument("--data_path",  action='store', default = data_file, dest = "data_path")
parser.add_argument("--save_folder",  action='store', default = './save_vae', dest = "save_folder")
parser.add_argument("--model",  default = "cvae", action='store', dest = "model")
parser.add_argument("--nz",  action='store', default = 3, type = int, dest = "nz")
parser.add_argument("--nb_epochs",  action='store', default = 1000, type = int, dest = "nb_epochs")
parser.add_argument("--batch_size",  action='store', default = 128, type = int, dest = "batch_size")
parser.add_argument("--lr",  action='store', default = 0.001, type = int, dest = "lr")

args = parser.parse_args()

data_path = args.data_path
save_folder = args.save_folder
model = args.model
nz = args.nz
batch_size = args.batch_size
nb_epochs = args.nb_epochs
lr = args.lr

###############################################################################
# Treating the data
# Imports the data
# Ensures that each trajectory has the same length (`padding`)
###############################################################################

data = pd.read_csv(data_path)
data.head()

lgv = (-4.28, 47.78)

padding = 24

traj_all = np.zeros((len(data.maree.unique()), 2, padding))

i = 0
for tt in data.maree.unique():
    traj = data[data.maree == tt].copy()
    traj = traj.loc[:,('x','y')].to_numpy()

    n = traj.shape[0]
    traj_all[i, :, 0:n+1] = traj.T
    i+=1

month_id = np.zeros((len(data.maree.unique())))
day_id = np.zeros((len(data.maree.unique())))
i = 0
for tt in data.maree.unique():
    month_id[i] = data[data.maree == tt].month.iloc[1]
    day_id[i] = datetime.weekday(datetime.fromisoformat(data[data.maree == tt].time.iloc[0]))
    i+=1

###############################################################################
# Only keeps weekdays trajectory
# Build covariates by projecting yday into a circle position
# Transforming longitude and latitude into standardized values centered around
# the fishing harbour position.
###############################################################################

traj_std = traj_all[day_id < 5,:,:].copy()

day = [datetime.strptime(data[data.maree == tt].time.iloc[1], '%Y-%m-%d').timetuple().tm_yday for tt in data.maree.unique()]
day = np.array(day)[day_id < 5]

day_x = np.cos(np.divide(day,365.25)*2*np.pi)
day_y = np.sin(np.divide(day,365.25)*2*np.pi)


# projection of coordinates in meters
R = 6377726 # earth radius
traj_std[:,0,:] = (np.pi * R * (traj_all[day_id < 5,0,:] + 4.28)/180)[:,:]
traj_std[:,1,:] = (np.pi * R * (traj_all[day_id < 5,1,:] - 47.78)/180)[:,:]

# standardization of (x,y) coordinates
scale = 10000 # empirical choice
traj_std = traj_std/scale

if model == 'cvae':
    traj_std = np.append(traj_std, np.array([day_x, day_y]).T[:,:,None], axis = 2)

###############################################################################
# Preparing the data and setting hyperparameters
# As the neural networks are used as a projection tools, the whole data serves
# as the train set
###############################################################################

class TrajDataSet(Dataset):
    def __init__(self, traj, transform = None):
        self.traj = traj
        self.transform = transform
    
    def __len__(self):
        return self.traj.shape[0]

    def __getitem__(self, idx):
        # select coordinates
        sample = self.traj[idx,:,:]
        if self.transform:
            sample = self.transform(sample)
        return sample

class ToTensor(object):
    """Convert ndarrays in sample to Tensors"""
    def __call__(self, sample):
        return torch.FloatTensor(sample)


def weights_init(m):
    """
    This function initializes the model weights randomly from a 
    Normal distribution. This follows the specification from the DCGAN paper.
    https://arxiv.org/pdf/1511.06434.pdf
    Source: https://pytorch.org/tutorials/beginner/dcgan_faces_tutorial.html
    """

    classname = m.__class__.__name__
    if classname.find('Conv') != -1:
        nn.init.normal_(m.weight.data, 0.0, 0.02)
    elif classname.find('BatchNorm') != -1:
        nn.init.normal_(m.weight.data, 1.0, 0.02)
        nn.init.constant_(m.bias.data, 0)

## hyperparameters
nb_traj = traj_std.shape[0]

## reduce size dataset
train_set = TrajDataSet(traj_std[:nb_traj,:,:], transform=ToTensor())
train_loader = DataLoader(train_set, batch_size=batch_size, 
    num_workers=0, shuffle=True, drop_last=False)

x = next(iter(train_loader))
x = x.cpu().numpy()



device = torch.device("cuda" if torch.cuda.is_available() else "cpu")


#torch.backends.cudnn.enabled = False

if model == 'cvae':
    E = cond_vae.cdepthCNNEncoder(nz =nz, nc = 2).to(device)
    G = cond_vae.cCNNGenerator(nz = nz, nc = 2).to(device)
if model == 'vae':
    E = cond_vae.CNNEncoder(nz =nz).to(device)
    G = cond_vae.CNNGenerator(nz = nz).to(device)

E = E.apply(weights_init)
G = G.apply(weights_init)

optimizer = 'adam'
lr = 0.001
optim_G = optim.Adam(G.parameters(), lr = lr, betas=(0.5, 0.999))
optim_E = optim.Adam(E.parameters(), lr = lr, betas=(0.5, 0.999))


niter = 0
traj_var = 1

monitor_dkl = []
monitor_loss = []
monitor_mse = []
monitor_traj_mse = []
best_loss = 1e12

for epoch in range(niter, nb_epochs+niter):
    print_dkl = []
    print_loss = []
    print_mse = []
    traj_mse = []

    for i, y in enumerate(train_loader):
        y = y.to(device)
        batch_size = y.shape[0]
        if model == 'cvae':
            x = y[:,:,0:padding]
            c = y[:,:,padding]
            c = c[:,:,None]
            Ex, mu, logVar = E(x, c)
            GEx = G(Ex, c)
            dkl = 0.5 * torch.sum(-1 - logVar + mu.pow(2) + logVar.exp())
            with torch.no_grad():
                Gmu = G(mu[:,:,None], c)
                l2 = 0.5 * F.mse_loss(x, Gmu, reduction='sum')
                traj_mse.append(l2.detach().item())
            loss_mse = 0.5 * F.mse_loss(x, GEx, reduction='sum')/(traj_var) 
        if model == 'vae':
            Ex, mu, logVar = E(y)
            GEx = G(Ex)           
            dkl = 0.5 * torch.sum(-1 - logVar + mu.pow(2) + logVar.exp())
            with torch.no_grad():
                Gmu = G(mu[:,:,None])
            with torch.no_grad():
                l2 = 0.5 * F.mse_loss(y, Gmu, reduction='sum')
                traj_mse.append(l2.detach().item())
            loss_mse = 0.5 * F.mse_loss(y, GEx, reduction='sum')/(traj_var) 

        loss = loss_mse + dkl
        print_dkl.append(dkl.detach().item())
        print_mse.append(loss_mse.detach().item())
        print_loss.append(loss.detach().item())


        optim_E.zero_grad()
        optim_G.zero_grad()

        loss.backward()

        optim_G.step()
        optim_E.step()
    monitor_dkl.append(np.sum(print_dkl)/nb_traj)
    monitor_loss.append(np.sum(print_loss)/nb_traj)
    monitor_mse.append(np.sum(print_mse)/nb_traj)
    monitor_traj_mse.append(np.sum(traj_mse)/nb_traj)
    if (epoch+1)%1 == 0:
        print('\nEpoch [{}/{}]: ----- [Tot/L2/KL/traj_L2]{:.4f}/{:.4f}/{:.4f}/{:.4f} -----  '
        .format(epoch + 1, nb_epochs, 
            np.sum(print_loss)/nb_traj, 
            np.sum(print_mse)/nb_traj, np.sum(print_dkl)/nb_traj,
            np.sum(traj_mse)/nb_traj))
    if np.sum(print_loss)/nb_traj < best_loss:
        best_loss = np.sum(print_loss)/nb_traj
        torch.save({
            'epoch': epoch,
            'E_state_dict': E.state_dict(),
            'G_state_dict': G.state_dict(),
            'optimizer_E_state_dict': optim_E.state_dict(),
            'optimizer_G_state_dict': optim_G.state_dict(),
            'loss': best_loss,
            'dkl': np.sum(print_dkl)/nb_traj,
            'rec_error': np.sum(traj_mse)/nb_traj
            }, f='%s/best_model_%d_%s_std_%s_epoch_%d.pth' % (save_folder, nz, model, optimizer, epoch+1))

    if (epoch+1)%500 == 0:
        torch.save(E.state_dict(), f='%s/Ed%d_%s_%s_epoch_%d.pth' % (save_folder, nz, model, optimizer, epoch+1))
        torch.save(G.state_dict(), f='%s/Gd%d_%s_%s_epoch_%d.pth' % (save_folder, nz, model, optimizer, epoch+1))

dict_monitor = {'loss' : monitor_loss, 'dkl' : monitor_dkl, 'mse' : monitor_mse, 'traj_mse' : monitor_traj_mse}
# open file for writing, "w" is writing
(pd.DataFrame.from_dict(data=dict_monitor)
   .to_csv("%s/monitor_%s_std_%d_%s_%d.csv" % (save_folder, nz, model, optimizer, nb_epochs), header=True))