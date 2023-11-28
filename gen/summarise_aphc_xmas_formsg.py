#!/usr/bin/env python3

# The MIT License
# Copyright (c) 2023 Adrian Tan <adrian_tan@nparks.gov.sg>
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
import random
import math


@click.command()
@click.option(
    "-w",
    "--working_dir",
    default=os.getcwd(),
    show_default=True,
    help="working directory",
)
def main(working_dir, phrase_file, latex_file):
    """
    generate bingo cards

    e.g. summarise_aphc_xmas_formsg.py
    """
    print("\t{0:<20} :   {1:<10}".format("working_dir", working_dir))
    print("\t{0:<20} :   {1:<10}".format("phrase_file", phrase_file))
    print("\t{0:<20} :   {1:<10}".format("latex_file", latex_file))

    # read phrase file
    phrases = []
    no_phrases = 0
    with open(phrase_file, "r") as file:
        for line in file:
            no_phrases += 1
            phrase = line.rstrip()
            phrases.append(phrase)

    print(f"{no_phrases} phrases read.")

    random.seed(10)
    stats = BingoStatistics()
    cards = []
    with open(latex_file, "w") as file:
        file.write(BingoMatrix.print_latex_header())
        for i in range(100):
            if i != 0:
                file.write("\\newpage\n")
            card = BingoMatrix(i + 1, 5, no_phrases)
            cards.append(card)
            stats.add_stats(card)
            file.write(card.print_latex(phrases))
        file.write(BingoMatrix.print_latex_footer())
    stats.report()

    # check overlap
    obs = []
    sum = 0
    sum2 = 0
    n = 1000000
    for i in range(n):
        x, y = random.sample(range(100), 2)
        # print(f'{x} {y}')
        overlap = 0
        for j in range(25):
            if cards[x].matrix[j] in cards[y].matrix:
                overlap += 1
        # print(f'\to={overlap}')
        sum += overlap
        obs.append(overlap)
    overlap_mean = sum / n
    resid2_sum = 0
    for i in range(n):
        resid2_sum += (obs[i] - overlap_mean) * (obs[i] - overlap_mean)
    std_dev = math.sqrt(resid2_sum / (n - 1))

    print(f"\n")
    print(f"overlap n   = {n}")
    print(f"expected overlap mean   = {25/52*25}")
    print(f"overlap mean            = {overlap_mean}")
    print(f"overlap std dev         = {std_dev:.2f}")
    print(f"\n")


class BingoStatistics(object):
    def __init__(self):
        self.no_cards = 0
        self.permutations = dict()
        self.phrases = dict()

    def add_card(self):
        self.no_cards += 1

    def add_phrase(self, phrase):
        if phrase in self.phrases.keys():
            self.phrases[phrase] += 1
        else:
            self.phrases[phrase] = 1

    def add_permutation(self, permutation):
        if permutation in self.permutations.keys():
            self.permutations[permutation] += 1
        else:
            self.permutations[permutation] = 1

    def add_stats(self, card):
        self.no_cards += 1
        diag1 = []
        diag2 = []
        for i in range(card.n):
            row = card.matrix[i * card.n : (i + 1) * card.n :]
            col = []
            for j in range(card.n):
                self.add_phrase(card.matrix[i * card.n + j])
                col.append(card.matrix[j * card.n + i])
            row.sort()
            col.sort()
            self.add_permutation(":".join(map(str, row)))
            self.add_permutation(":".join(map(str, col)))
            diag1.append(card.matrix[i * card.n + i])
            diag2.append(card.matrix[(card.n - i - 1) * card.n + i])
        diag1.sort()
        diag2.sort()
        self.add_permutation(":".join(map(str, diag1)))
        self.add_permutation(":".join(map(str, diag2)))

    def report(self):
        print("phrases")
        k = list(self.phrases.keys())
        k.sort()
        sum = 0
        sumX2 = 0
        n = len(k)
        for key in k:
            sum += self.phrases[key]
            sumX2 += self.phrases[key] * self.phrases[key]

        print(f"no of occurrence of phrases        : {sum}")
        print(f"no of unique phrases               : {n}")
        print(f"mean occurrence of phrases         : {sum/n:.2f}")
        print(
            f"std dev of occurrence of phrases   : {math.sqrt(sumX2/n-(sum/n*sum/n)):.2f}"
        )
        print(f"expected mean                      : {2500/52:.4f}")
        print(f"no of unique permutations          : {len(self.permutations.keys())}")
        print(f"no of rows, cols and diags         : {self.no_cards*12}")

    def print(self):
        print("phrases")
        print(self.phrases)
        print("permutations")
        print(self.permutations)


class BingoMatrix(object):
    def __init__(self):
        self.matrix = []

    def __init__(self, id, n, s):
        self.id = id
        self.n = n
        self.s = s
        self.matrix = random.sample(range(s), n * n)

    def print(self):
        for i in range(self.n):
            for j in range(self.n):
                print(f"  {self.matrix[i*self.n + j]:2d}", end="")
            print("\n")
        print("\n")

    @staticmethod
    def print_latex_header():
        text = "\documentclass[20pt,A4]{article}\n"
        text += "\\usepackage{array}\n"
        text += (
            "\\usepackage[centering, margin={1in, 0.8in}, includeheadfoot]{geometry}\n"
        )
        text += "\\usepackage{graphicx}\n"
        text += "\\usepackage[table]{xcolor}\n"
        text += "\\usepackage{mathptmx}\n"
        text += "\\usepackage{anyfontsize}\n"
        text += "\\usepackage{t1enc}\n"
        text += "\\usepackage{background}\n"
        text += "\\usepackage[absolute,overlay]{textpos}\n"
        text += "\n"
        text += "\setlength{\headheight}{0in} \n"
        text += "\setlength{\headsep}{0pt}\n"
        text += "\setlength{\\topmargin}{0pt}\n"
        text += "\setlength{\\footskip}{0pt}\n"
        text += "\setlength{\marginparsep}{11pt}\n"
        text += "\pagestyle{empty}\n"
        text += "\n"
        text += "\\backgroundsetup{contents=\includegraphics{cny3.jpg}, hshift=0cm, vshift=1cm, angle=0, scale=0.38, opacity=0.85}\n"
        text += "\\newcolumntype{L}[1]{>{\\raggedright\let\\newline\\\\\\arraybackslash\hspace{0pt}}m{#1}}\n"
        text += "\\newcolumntype{C}[1]{>{\centering\let\\newline\\\\\\arraybackslash\hspace{0pt}}m{#1}}\n"
        text += "\\newcolumntype{R}[1]{>{\\raggedleft\let\\newline\\\\\\arraybackslash\hspace{0pt}}m{#1}}\n"
        text += "\n"
        text += "\definecolor{darkred}{rgb}{1,0.35,0.25}\n"
        text += "\n"
        text += "\\begin{document}\n"
        return text

    @staticmethod
    def print_latex_footer():
        text = "\end{document}\n"
        return text

    def print_latex(self, phrases):
        text = ""
        text += "\\begin{textblock*}{3cm}(15.5cm,4cm)\n"
        text += "     \Huge \\textcolor{darkred}{" + f"{self.id:2d}" + "}\n"
        text += "\end{textblock*}\n"
        text += "\n"
        text += "\hphantom\n"
        text += "\\newline\n"
        text += "\\newline\n"
        text += "\\newline\n"
        text += "\\newline\n"
        text += "\\newline\n"
        text += "\\newline\n"
        text += "\\newline\n"
        text += "\\newline\n"
        text += "\n"
        text += "\\begin{center}\n"
        text += "\large \\bfseries\n"
        text += "\\begin{tabular}{| C{2.8cm}| C{2.8cm}| C{2.8cm}| C{2.8cm}| C{2.8cm}|} \hline  \n"
        for i in range(self.n):
            text += f"{phrases[self.matrix[i*self.n]]} & "
            text += f"{phrases[self.matrix[i*self.n+1]]} & "
            text += f"{phrases[self.matrix[i*self.n+2]]} & "
            text += f"{phrases[self.matrix[i*self.n+3]]} & "
            text += f"{phrases[self.matrix[i*self.n+4]]} \\\\[2cm]\n"
            text += "\hline\n"
        text += "\end{tabular}\n"
        text += "\end{center}\n"
        text += "\n"
        text += "\centerline{\\textcolor{darkred}{brought to you by the APHC Bonding Committee}}\n"
        return text


if __name__ == "__main__":
    main()
