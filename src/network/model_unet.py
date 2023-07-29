# code from https://gitee.com/LynnHg/pytorch-medical-image-segmentation

import sys
import torch
import torch.nn as nn
import torch.nn.functional as F

activation = nn.ReLU

def initialize_weights(*models, a=0):
    for model in models:
        for module in model.modules():
            if isinstance(module, nn.Conv2d) or isinstance(module, nn.Linear) or isinstance(module, nn.ConvTranspose2d):
                nn.init.kaiming_normal_(module.weight, a=a)
                if module.bias is not None:
                    nn.init.constant_(module.bias, 0)
            elif isinstance(module, nn.BatchNorm2d):
                module.weight.data.fill_(1)
                module.bias.data.zero_()


class EncoderBlock(nn.Module):
    def __init__(self, ch_in, ch_out, depth=2, use_res=False):
        super(EncoderBlock, self).__init__()

        self.use_res = use_res

        self.conv = [nn.Sequential(
            nn.Conv2d(ch_in, ch_out, kernel_size=3, padding=1, bias=False),
            nn.BatchNorm2d(ch_out),
            activation(),
        )]

        for i in range(1, depth):
            self.conv.append(nn.Sequential(nn.Conv2d(ch_out, ch_out, kernel_size=3, padding=1, bias=False),
                                           nn.BatchNorm2d(ch_out),
                                           nn.Sequential() if use_res and i == depth-1 else activation()
                                           ))
        self.conv = nn.Sequential(*self.conv)
        if use_res:
            self.conv1x1 = nn.Conv2d(ch_in, ch_out, 1)

    def forward(self, x):
        if self.use_res:
            residual = self.conv1x1(x)

        x = self.conv(x)

        if self.use_res:
            x += residual
            x = F.relu(x)

        return x


class DecoderBlock(nn.Module):
    """
    Interpolate
    """

    def __init__(self, ch_in, ch_out, use_deconv=False):
        super(DecoderBlock, self).__init__()
        if use_deconv:
            self.up = nn.Sequential(
                nn.ConvTranspose2d(ch_in, ch_out, kernel_size=2, stride=2)
            )
        else:
            self.up = nn.Sequential(
                nn.Upsample(scale_factor=2, mode='bilinear', align_corners=False),
                nn.Conv2d(ch_in, ch_out, kernel_size=3, padding=1, bias=False),
                nn.BatchNorm2d(ch_out),
                activation()
            )

    def forward(self, x):
        return self.up(x)


class Unet(nn.Module):
    def __init__(self, img_ch=3, num_classes=5, depth=2):
        super(Unet, self).__init__()

        chs = [64, 128, 256, 512, 512]

        self.maxpool = nn.MaxPool2d(kernel_size=2, stride=2)

        self.enc1 = EncoderBlock(img_ch, chs[0], depth=depth)
        self.enc2 = EncoderBlock(chs[0], chs[1], depth=depth)
        self.enc3 = EncoderBlock(chs[1], chs[2], depth=depth)
        self.enc4 = EncoderBlock(chs[2], chs[3], depth=depth)
        self.enc5 = EncoderBlock(chs[3], chs[4], depth=depth)

        self.dec4 = DecoderBlock(chs[4], chs[3])
        self.decconv4 = EncoderBlock(chs[3] * 2, chs[3])

        self.dec3 = DecoderBlock(chs[3], chs[2])
        self.decconv3 = EncoderBlock(chs[2] * 2, chs[2])

        self.dec2 = DecoderBlock(chs[2], chs[1])
        self.decconv2 = EncoderBlock(chs[1] * 2, chs[1])

        self.dec1 = DecoderBlock(chs[1], chs[0])
        self.decconv1 = EncoderBlock(chs[0] * 2, chs[0])

        self.conv_1x1 = nn.Conv2d(chs[0], num_classes, 1, bias=False)

        initialize_weights(self)

    def forward(self, x):
        # encoding path
        x1 = self.enc1(x)
        x2 = self.maxpool(x1)

        x2 = self.enc2(x2)
        x3 = self.maxpool(x2)

        x3 = self.enc3(x3)
        x4 = self.maxpool(x3)

        x4 = self.enc4(x4)
        x5 = self.maxpool(x4)

        x5 = self.enc5(x5)

        # decoding + concat path
        d4 = self.dec4(x5)
        d4 = torch.cat((x4, d4), dim=1)
        d4 = self.decconv4(d4)

        d3 = self.dec3(d4)
        d3 = torch.cat((x3, d3), dim=1)
        d3 = self.decconv3(d3)

        d2 = self.dec2(d3)
        d2 = torch.cat((x2, d2), dim=1)
        d2 = self.decconv2(d2)

        d1 = self.dec1(d2)
        d1 = torch.cat((x1, d1), dim=1)
        d1 = self.decconv1(d1)

        d1 = self.conv_1x1(d1)

        return d1
