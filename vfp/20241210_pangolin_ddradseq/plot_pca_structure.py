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
@click.option("-o", "--output_dir", required=False, default=os.path.join(os.getcwd(), f"plot-{uuid.uuid4()}"), help="output directory", type=str)
@click.option("-z", "--output_pdf_file", required=False, help="output pdf file", type=str)
@click.option("-x", "--x_axis_header", required=False, help="x axis", type=str)
@click.option("-y", "--y_axis_header", required=True, help="y axis", type=str)
@click.option("-t", "--title", required=True, help="title", type=str)
@click.option("-p", "--point", default=1, help="Point character - http://www.phaget4.org/R/plot.html", type=int)
@click.option("-q", "--cex", default=1, help="Scale Point character [0,1]", type=float)
@click.option("-c", "--colour", default="black", help="Point colour - black, red, blue, green", type=str)
@click.option("--ylim", help="Plot y axis on [0,1]", type=str)
@click.option("-k", "--keep", default=False, is_flag=True, help="To keep output directory or not")
def main(input_txt_file, output_dir, output_pdf_file, x_axis_header, y_axis_header, title, point, cex, colour, ylim, keep):
    """
    qwikplot [options] file

    -o   output directory
    -z   output pdf file name
    -x   x header name
    -y   y header name
    -s   sort y axis values
    -t   Title
    -u   plot [0,1] on both axes
    -p   Point Character: 0-25 (Default 1)
         http://www.phaget4.org/R/plot.html
    -c   Colour : black, red, blue, green (Default black)
    -q   scale points [0,1]

    e.g. qwikplot sample.txt -y call_rate -t "call rate"
         qwikplot sample.txt -x PC1 -y PC2 -t "PC1 vs PC2" -c blue -p 4
    """
    print("\t{0:<20} :   {1:<10}".format("file", input_txt_file))
    print("\t{0:<20} :   {1:<10}".format("output dir", output_dir))
    print("\t{0:<20} :   {1:<10}".format("output pdf", output_pdf_file))
    if x_axis_header is not None:
        print("\t{0:<20} :   {1:<10}".format("x axis header", x_axis_header))
    print("\t{0:<20} :   {1:<10}".format("y axis header", y_axis_header))
    print("\t{0:<20} :   {1:<10}".format("title", title))
    print("\t{0:<20} :   {1:<10}".format("point", point))
    print("\t{0:<20} :   {1:<10}".format("colour", colour))
    print("\t{0:<20} :   {1:<10}".format("keep", keep))

    #programs
    rscript = "/usr/bin/Rscript"

    try:
        os.makedirs(output_dir, exist_ok=True)
    except OSError as error:
        print(f"{error.filename} cannot be created")

    input_txt_file = os.path.abspath(input_txt_file)

    if output_pdf_file is None:
        if x_axis_header is not None:
            output_pdf_file = f"{output_dir}/{x_axis_header}-vs-{y_axis_header}.pdf"
        else:
            output_pdf_file = f"{output_dir}/{y_axis_header}.pdf"
    else:
        if os.path.isabs(output_pdf_file) is False:
            output_pdf_file = os.path.join(os.getcwd(), output_pdf_file)

    #generate R script
    r_script_file =f"{output_dir}/plot.r"
    if ylim is not None:
        ylim = f"ylim={ylim},"
    else:
        ylim = ""
    if x_axis_header is None:
        with open(r_script_file, "w") as f:
            r_script = f"""
            data = read.table("{input_txt_file}", header = TRUE, comment.char="")
            colnames(data) <- gsub("^X.", "", colnames(data))
            pdf("{output_pdf_file}")
            plot(sort(data${y_axis_header}), type = "p", cex={cex}, pch={point}, col = "{colour}", {ylim} xlab = "index", ylab = "{y_axis_header}", main = "{title}")
            dev.off()"""
            f.write(r_script)
    else:
        with open(r_script_file, "w") as f:
            r_script = f"""
            data = read.table("{input_txt_file}", header = TRUE, comment.char="")
            colnames(data) <- gsub("^X\\.", "", colnames(data))
            pdf("{output_pdf_file}")
            plot(data${x_axis_header}, data${y_axis_header}, type = "p", cex={cex}, pch = {point}, col = "{colour}", {ylim} xlab = "{x_axis_header}", ylab = "{y_axis_header}", main = "{title}")
            dev.off()"""
            f.write(r_script)

    #run R script
    subprocess.run(f"{rscript} {r_script_file}", shell=True, check=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)

    if not keep:
        subprocess.run(f"rm -fr {output_dir}", shell=True, check=True)

class MiniPipeManager(object):
    def __init__(self, log_file):
        self.log_file = log_file
        self.log_msg = []

    def run(self, cmd, tgt, desc):
        try:
            if os.path.exists(tgt):
                self.log(f"{desc} -  already executed")
                self.log(cmd)
                return
            else:
                self.log(f"{desc}")
                self.log(cmd)
                subprocess.run(cmd, shell=True, check=True)
                subprocess.run(f"touch {tgt}", shell=True, check=True)
        except subprocess.CalledProcessError as e:
            self.log(f" - failed")
            exit(1)

    def log(self, msg):
        print(msg)
        self.log_msg.append(msg)

    def print_log(self):
        self.log(f"\nlogs written to {self.log_file}")
        with open(self.log_file, "w") as f:
            f.write("\n".join(self.log_msg))


if __name__ == "__main__":
    main() # type: ignore