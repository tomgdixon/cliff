#!/usr/bin/env python3
"""
Generate independent ground-truth Fast Fourier Transform (FFT) test vectors
using standard NumPy FFT (numpy.fft.fft, numpy.fft.fft2) to cross-validate
against Haskell's Clifford Fast Fourier Transform (fcft1D, fcft2D).

In Cl(2,0,0), even multivectors a + b*e12 are isomorphic to complex numbers a + b*i,
and the Clifford bivector kernel exp(-e12 * theta) matches the complex exponential exp(-i * theta).
"""

import numpy as np

def generate_haskell_code():
    np.random.seed(42)
    
    # 1. 1D FFT test cases: sizes 4, 8, 16, 64
    sizes_1d = [4, 8, 16, 64]
    cases_1d = []
    
    for n in sizes_1d:
        for trial in range(3):
            # Generate random complex signal
            real = np.random.uniform(-5.0, 5.0, n)
            imag = np.random.uniform(-5.0, 5.0, n)
            signal = real + 1j * imag
            
            # NumPy Forward FFT
            fft_out = np.fft.fft(signal)
            
            # Pack as (real_in, imag_in, real_out, imag_out)
            r_in = [float(x) for x in real]
            i_in = [float(x) for x in imag]
            r_out = [float(x) for x in fft_out.real]
            i_out = [float(x) for x in fft_out.imag]
            
            cases_1d.append((n, r_in, i_in, r_out, i_out))
            
    # 2. 2D FFT test cases: shapes (4, 4) and (8, 8)
    shapes_2d = [(4, 4), (8, 8)]
    cases_2d = []
    
    for (h, w) in shapes_2d:
        for trial in range(2):
            real = np.random.uniform(-5.0, 5.0, (h, w))
            imag = np.random.uniform(-5.0, 5.0, (h, w))
            grid = real + 1j * imag
            
            # NumPy 2D FFT
            fft2_out = np.fft.fft2(grid)
            
            r_in = [float(x) for x in real.flatten()]
            i_in = [float(x) for x in imag.flatten()]
            r_out = [float(x) for x in fft2_out.real.flatten()]
            i_out = [float(x) for x in fft2_out.imag.flatten()]
            
            cases_2d.append((h, w, r_in, i_in, r_out, i_out))
            
    # Format as Haskell data
    lines = []
    lines.append("-- Auto-generated FFT oracle test cases from NumPy 1D and 2D FFT")
    lines.append("fft1DCases :: [(Int, [Double], [Double], [Double], [Double])]")
    lines.append("fft1DCases =")
    lines.append("  [")
    for idx, (n, ri, ii, ro, io) in enumerate(cases_1d):
        sep = "," if idx < len(cases_1d) - 1 else ""
        lines.append(f"    ( {n}")
        lines.append(f"    , {ri}")
        lines.append(f"    , {ii}")
        lines.append(f"    , {ro}")
        lines.append(f"    , {io}")
        lines.append(f"    ){sep}")
    lines.append("  ]")
    lines.append("")
    lines.append("fft2DCases :: [(Int, Int, [Double], [Double], [Double], [Double])]")
    lines.append("fft2DCases =")
    lines.append("  [")
    for idx, (h, w, ri, ii, ro, io) in enumerate(cases_2d):
        sep = "," if idx < len(cases_2d) - 1 else ""
        lines.append(f"    ( {h}, {w}")
        lines.append(f"    , {ri}")
        lines.append(f"    , {ii}")
        lines.append(f"    , {ro}")
        lines.append(f"    , {io}")
        lines.append(f"    ){sep}")
    lines.append("  ]")
    
    return "\n".join(lines)

if __name__ == "__main__":
    code = generate_haskell_code()
    with open("test/Test/FFTOracleData.hs", "w") as f:
        f.write("module Test.FFTOracleData (fft1DCases, fft2DCases) where\n\n" + code + "\n")
    print("Successfully generated test/Test/FFTOracleData.hs")
