from pathlib import Path
import random

ROOT = Path(__file__).resolve().parents[1] / "sample_data" / "genomes"

def wrap(seq, n=70): return "\n".join(seq[i:i+n] for i in range(0,len(seq),n))+"\n"
def motif(i, length=41):
    rng = random.Random(1000 + i)
    return "ATG" + "".join(rng.choice("ACGT") for _ in range(length-6)) + "TAA"

def main():
    a_dir,b_dir=ROOT/"RefA",ROOT/"RefB"; a_dir.mkdir(parents=True,exist_ok=True); b_dir.mkdir(parents=True,exist_ok=True)
    spacer=lambda n,seed: "".join("ACGT"[(i*5+seed*3+i//7)%4] for i in range(n))
    genes=[motif(i) for i in range(1,9)]
    b=spacer(99,1); b_pos=[]
    for i,g in enumerate(genes):
        b_pos.append(len(b)+1); b+=g+spacer(54,i+3)
    b+=spacer(150,8)
    # Same genes in forward order at a shifted/larger target interval.
    a=spacer(169,9); a_pos=[]
    for i,g in enumerate(genes):
        a_pos.append(len(a)+1); a+=g+spacer(69,i+12)
    a+=spacer(170,19)
    (a_dir/"refA.fa").write_text(">Chr1 canonical artificial genome\n"+wrap(a),encoding="ascii")
    (b_dir/"refB.fa").write_text(">Chr1 source artificial genome\n"+wrap(b),encoding="ascii")
    for directory,positions in ((a_dir,a_pos),(b_dir,b_pos)):
        rows=["##gff-version 3"]
        for i,pos in enumerate(positions,1):
            end=pos+len(genes[i-1])-1; gid=f"gene{i}"
            rows += [f"Chr1\tQTLift\tgene\t{pos}\t{end}\t.\t+\t.\tID={gid};Name={gid}",f"Chr1\tQTLift\tmRNA\t{pos}\t{end}\t.\t+\t.\tID={gid}.t1;Parent={gid}",f"Chr1\tQTLift\tCDS\t{pos}\t{end}\t.\t+\t0\tID={gid}.cds;Parent={gid}.t1"]
        (directory/("refA.gff3" if directory==a_dir else "refB.gff3")).write_text("\n".join(rows)+"\n",encoding="ascii")
    print(ROOT)
if __name__=="__main__": main()
