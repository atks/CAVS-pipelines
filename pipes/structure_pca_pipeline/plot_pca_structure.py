#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2025 Adrian Tan <adrian_tan@nparks.gov.sg>
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the 'Software'), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED 'AS IS', WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import os
import click
import subprocess
import uuid

@click.command()
@click.argument("input_txt_file", type=click.Path(exists=True))
@click.option("-o", "--output_r_file", required=True, help="output R file", type=str)
@click.option("-z", "--output_pdf_file", required=True, help="output pdf file", type=str)
def main(input_txt_file, output_r_file, output_pdf_file):
    """
    plot_pca_structure.py [options] <input_txt_file>

    -o   output directory
    -z   output pdf file name

    e.g. plot_pca_structure.py sample.txt -o gisplot -z "k2_pangolin.pdf"
    """
    print("\t{0:<20} :   {1:<10}".format("file", input_txt_file))
    print("\t{0:<20} :   {1:<10}".format("output R file", output_r_file))
    print("\t{0:<20} :   {1:<10}".format("output pdf", output_pdf_file))

    #programs
    rscript = "/usr/bin/Rscript"

    input_txt_file = os.path.abspath(input_txt_file)

    #generate R script
    r_script_file = output_r_file

    with open(r_script_file, "w") as f:
        r_script = f"""
        library(ggplot2)
        library(scatterpie)
        data <- read.table("{input_txt_file}", header = T)
        K = ncol(data) - 21
        data$region <- factor(1:nrow(data))
        COLS <- paste0("C", 1:K)
        pdf("{output_pdf_file}")
        ggplot() +
            geom_scatterpie(aes(x = PC1, y = PC2, group = region), data = data, cols = COLS[1:K]) +
            coord_equal()
        dev.off()"""
        f.write(r_script)

    #run R script
    subprocess.run(f"{rscript} {r_script_file}", shell=True, check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

if __name__ == "__main__":
    main() # type: ignore