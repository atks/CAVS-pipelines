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
@click.option("-o", "--output_dir", required=True, default=os.path.join(os.getcwd(), f"plot-{uuid.uuid4()}"), help="output directory", type=str)
@click.option("-z", "--output_pdf_file", required=True, help="output pdf file", type=str)
def main(input_txt_file, output_dir, output_pdf_file):
    """
    plot_gis_structure.py [options] <input_txt_file>

    -o   output directory
    -z   output pdf file name

    e.g. plot_gis_structure.py sample.txt -o gisplot -z "k2_pangolin.pdf"
    """
    print("\t{0:<20} :   {1:<10}".format("file", input_txt_file))
    print("\t{0:<20} :   {1:<10}".format("output dir", output_dir))
    print("\t{0:<20} :   {1:<10}".format("output pdf", output_pdf_file))

    #programs
    rscript = "/usr/bin/Rscript"

    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    input_txt_file = os.path.abspath(input_txt_file)

    #generate R script
    r_script_file =f"{output_dir}/plot.r"

    with open(r_script_file, "w") as f:
        r_script = f"""
        library(ggplot2)
        library(ggmap)
        library(scatterpie)
        data <- read.table("{input_txt_file}", header = T)
        K = ncol(data) - 4
        data$region <- factor(1:nrow(data))
        head(data)
        register_google(key = "AIzaSyDxpiYURcyH35qE2PprCh3YVZAqy7yPxtw")
        singapore_city <- geocode("Singapore")
        singapore_map <- get_map(singapore_city, zoom = 11, maptype = "satellite")
        COLS <- paste0("C", 1:K)
        pdf("{output_pdf_file}")
        ggmap(singapore_map) +
            geom_scatterpie(aes(x = longitude, y = latitude, group = region), data = data, cols = COLS[1:K]) +
            coord_equal()
            dev.off()"""
        f.write(r_script)

    #run R script
    subprocess.run(f"{rscript} {r_script_file}", shell=True, check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

if __name__ == "__main__":
    main() # type: ignore