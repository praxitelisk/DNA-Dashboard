import streamlit as st
import io
import neatbio.sequtils as utils
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
import numpy as np

st.set_option('deprecation.showfileUploaderEncoding', False)

def main():
    st.title('Simple Bioinformatics App')
    menu = ["DNA Sequence", "Dot Plot"]
    choice = st.sidebar.selectbox("Select Activity", menu)

    if choice == "DNA Sequence":
        st.subheader("DNA Sequence Analysis")
        seq_file = st.file_uploader("Upload FASTA File", type=["fasta", "fa"])

        if seq_file is not None:
            # Convert binary file to text
            stringio = io.StringIO(seq_file.getvalue().decode("utf-8"))
            dna_record = SeqIO.read(stringio, "fasta")
            dna_seq = dna_record.seq
            desc = dna_record.description

            details = st.radio("Details", ("Description", "Sequence"))
            if details == "Description":
                st.write(desc)
            else:
                st.write(dna_seq)

            # Nucleotide Frequencies
            st.subheader("Nucleotide Frequencies")
            dna_freq = Counter(dna_seq)
            st.write(dna_freq)

            # Color pickers for nucleotides
            adenine_color = st.color_picker("Adenine (A) Color", "#FF0000")
            guanine_color = st.color_picker("Guanine (G) Color", "#00FF00")
            thymine_color = st.color_picker("Thymine (T) Color", "#0000FF")
            cytosine_color = st.color_picker("Cytosine (C) Color", "#FFFF00")

            if st.button("Plot Nucleotide Frequencies"):
                fig, ax = plt.subplots()
                bars = ax.bar(dna_freq.keys(), dna_freq.values())
                bars[0].set_color(adenine_color)  # A
                bars[1].set_color(thymine_color)  # T
                bars[2].set_color(guanine_color)  # G
                bars[3].set_color(cytosine_color) # C
                st.pyplot(fig)

            # DNA Composition
            st.subheader("DNA Composition")
            gc_score = utils.gc_content(str(dna_seq))
            at_score = utils.at_content(str(dna_seq))
            st.json({"GC Content": gc_score, "AT Content": at_score})

            # Nucleotide Count
            nt_count = st.text_input("Enter Nucleotide Here", "").upper()
            if nt_count:
                st.write(f"Number of {nt_count} nucleotides: {str(dna_seq).count(nt_count)}")

            # Protein Synthesis
            p1 = dna_seq.translate()
            aa_freq = Counter(str(p1))

            st.subheader("Protein Synthesis")
            if st.checkbox("Transcription"):
                st.write(dna_seq.transcribe())
            if st.checkbox("Translation"):
                st.write(dna_seq.translate())
            if st.checkbox("Complement"):
                st.write(dna_seq.complement())
            if st.checkbox("Amino Acid Frequency"):
                st.write(aa_freq)
            if st.checkbox("Amino Acid Frequency Plot"):
                fig, ax = plt.subplots()
                ax.bar(aa_freq.keys(), aa_freq.values())
                st.pyplot(fig)
            if st.checkbox("Full Amino Acid Name"):
                aa_name = str(p1).replace("*", "")
                st.write(aa_name)
                st.write("--------------------------")
                st.write(utils.convert_1to3(aa_name))

    elif choice == "Dot Plot":
        st.subheader("Generate Dot Plot For Two Sequences")
        seq_file1 = st.file_uploader("Upload 1st FASTA File", type=["fasta", "fa"])
        seq_file2 = st.file_uploader("Upload 2nd FASTA File", type=["fasta", "fa"])

        if seq_file1 and seq_file2:
            # Convert binary files to text
            stringio1 = io.StringIO(seq_file1.getvalue().decode("utf-8"))
            stringio2 = io.StringIO(seq_file2.getvalue().decode("utf-8"))
            dna_record1 = SeqIO.read(stringio1, "fasta")
            dna_record2 = SeqIO.read(stringio2, "fasta")

            dna_seq1 = dna_record1.seq
            dna_seq2 = dna_record2.seq

            desc1 = dna_record1.description
            desc2 = dna_record2.description

            details = st.radio("Details", ("Description", "Sequence"))
            if details == "Description":
                st.write(desc1)
                st.write("----------")
                st.write(desc2)
            else:
                st.write(dna_seq1)
                st.write("----------")
                st.write(dna_seq2)

            custom_limit = st.number_input("Select max number of Nucleotides", 10, 200, 25)
            if st.button("Generate Dot Plot"):
                st.write(f"Comparing the first {custom_limit} nucleotides of two sequences")
                fig = dotplotx(dna_seq1[:custom_limit], dna_seq2[:custom_limit])
                st.pyplot(fig)

def delta(x, y):
    return 0 if x == y else 1

def makeMatrix(seq1, seq2, k):
    n, m = len(seq1), len(seq2)
    return [[sum(delta(x, y) for x, y in zip(seq1[i:i+k], seq2[j:j+k])) for j in range(m-k+1)] for i in range(n-k+1)]

def dotplotx(seq1, seq2):
    matrix = np.array(makeMatrix(seq1, seq2, 1))
    fig, ax = plt.subplots()
    ax.imshow(matrix, cmap="Greys", interpolation="nearest")
    ax.set_xticks(np.arange(len(seq2)))
    ax.set_xticklabels(list(seq2))
    ax.set_yticks(np.arange(len(seq1)))
    ax.set_yticklabels(list(seq1))
    return fig

if __name__ == '__main__':
    main()
