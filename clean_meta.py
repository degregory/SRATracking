#!/bin/env python3

import os
import sys

lines = []
with open("SRA_meta.tsv", "r") as in_fh:
    lines = in_fh.readlines()

print(len(lines))
lines = list(set(lines))

print(len(lines))

with open("SRA_meta.tsv", "w") as out_fh:
    out_fh.writelines(lines)