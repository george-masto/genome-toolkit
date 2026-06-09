#!/usr/bin/env python3
"""Generate the README banner — an original DNA double-helix graphic.

Reproducible, dependency-light (matplotlib + numpy), and commit-safe (it's our
own generated art, no personal data). Run:  python3 assets/banner/make_banner.py
Outputs: assets/banner.png
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import patheffects

W, H = 1600, 460
fig = plt.figure(figsize=(W / 100, H / 100), dpi=100)
ax = fig.add_axes([0, 0, 1, 1])
ax.axis("off")

# --- diagonal gradient background (deep navy -> teal) ---
xx, yy = np.meshgrid(np.linspace(0, 1, W), np.linspace(0, 1, H))
t = 0.6 * xx + 0.4 * yy
c1 = np.array([0x0d, 0x11, 0x26]) / 255.0   # deep navy
c2 = np.array([0x10, 0x47, 0x56]) / 255.0   # teal
bg = c1[None, None, :] + (c2 - c1)[None, None, :] * t[:, :, None]
ax.imshow(bg, extent=[0, W, 0, H], origin="lower", aspect="auto", zorder=0)

# --- DNA double helix ---
xs = np.linspace(0, 4 * np.pi, 1200)
amp, mid = H * 0.17, H * 0.5
px = xs * (W / (4 * np.pi))
s1 = mid + amp * np.sin(xs)
s2 = mid + amp * np.sin(xs + np.pi)
ax.plot(px, s1, color="#5ad1c8", lw=3, alpha=0.45, zorder=1)
ax.plot(px, s2, color="#7aa2ff", lw=3, alpha=0.45, zorder=1)
cmap = plt.cm.cool
for i in range(0, len(xs), 26):
    ax.plot([px[i], px[i]], [s1[i], s2[i]],
            color=cmap((np.sin(xs[i]) + 1) / 2), lw=2, alpha=0.35, zorder=1)

# --- title + tagline ---
title = ax.text(W * 0.5, H * 0.60, "genome-toolkit", color="white",
                fontsize=66, fontweight="bold", ha="center", va="center", zorder=3)
title.set_path_effects([patheffects.withStroke(linewidth=5, foreground=(0, 0, 0, 0.55))])
sub = ax.text(W * 0.5, H * 0.33,
              "your 30× whole genome — from raw reads to insight, privately",
              color="#d4ecff", fontsize=23, ha="center", va="center", zorder=3)
sub.set_path_effects([patheffects.withStroke(linewidth=3, foreground=(0, 0, 0, 0.45))])

ax.set_xlim(0, W)
ax.set_ylim(0, H)
out = os.path.join(os.path.dirname(__file__), "..", "banner.png")
fig.savefig(out, dpi=100)
print("wrote", os.path.normpath(out))
