import argparse

import numpy as np
import torch
import SimpleITK as sitk
import image_viewing
from GUI.viewer_utils import ImageMask

parser = argparse.ArgumentParser(description='3D Slice Image Viewer.')
parser.add_argument('img', type=str, help='path to input image')
parser.add_argument('-l', type=str, help='path to label image')
parser.add_argument('-p', type=str, help='path to landmarks (torch tensor object or csv file)')
parser.add_argument('--pformat', type=str, choices=['xyz', 'zyx', 'xzy', 'zxy', 'yxz', 'yzx'], default='xyz',
                    help='format of landmarks, i.e. the interpretation of coordinates')
parser.add_argument('--pvalues', type=str, choices=['unit', 'coord'], default='coord',
                    help='either points are converted to unit sphere (pytorch representation) or they are coordinates corresponding to the image')
args = parser.parse_args()

img = sitk.ReadImage(args.img)
if args.l:
    lbl = sitk.ReadImage(args.l)
    lbl = ImageMask(lbl)
else:
    lbl = None

if args.p:
    if args.p.endswith('.csv'):
        pts = torch.from_numpy(np.loadtxt(args.p, delimiter=','))
    elif args.p.endswith('.pt') or args.p.endswith('.pth'):
        pts = torch.load(args.p, map_location=torch.device('cpu'))
    else:
        raise ValueError('Unknown file format for landmarks file: {}'.format(args.p))
    
    # ensure shape of points (N x 3)
    if pts.shape[0] == 3:
        pts = pts.transpose(0, 1)
    
    # ensure xyz format
    if args.pformat != 'xyz':
        permutation = (args.pformat.find('x'), args.pformat.find('y'), args.pformat.find('z'))
        pts = torch.stack([pts[:, permutation[0]], pts[:, permutation[1]], pts[:, permutation[2]]], dim=1)

    # ensure coordinate values
    if args.pvalues == 'unit':
        pts = (pts + 1) / 2 * torch.tensor(img.GetSize())

    pts = pts.tolist()
else:
    pts = None

image_viewing._show_image(img, window_title=args.img, blocking=True, mask=lbl, return_markers=False, points=pts)
