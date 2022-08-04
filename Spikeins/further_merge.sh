outdir="WT/1-seq"
mkdir -p  $outdir
cat  A1_S13/A1_S13.spikein-1.seq   A2_S14/A2_S14.spikein-1.seq  >   $outdir/spikein-1.seq
cat  A1_S13/A1_S13.spikein-2.seq   A2_S14/A2_S14.spikein-2.seq  >   $outdir/spikein-2.seq
cat  A1_S13/A1_S13.spikein-3.seq   A2_S14/A2_S14.spikein-3.seq  >   $outdir/spikein-3.seq
cat  A1_S13/A1_S13.spikein-4.seq   A2_S14/A2_S14.spikein-4.seq  >   $outdir/spikein-4.seq
cat  A1_S13/A1_S13.spikein-5.seq   A2_S14/A2_S14.spikein-5.seq  >   $outdir/spikein-5.seq


outdir="INPUT/1-seq"
mkdir -p  $outdir
cat  Ar_S31/Ar_S31.spikein-1.seq      >   $outdir/spikein-1.seq
cat  Ar_S31/Ar_S31.spikein-2.seq      >   $outdir/spikein-2.seq
cat  Ar_S31/Ar_S31.spikein-3.seq      >   $outdir/spikein-3.seq
cat  Ar_S31/Ar_S31.spikein-4.seq      >   $outdir/spikein-4.seq
cat  Ar_S31/Ar_S31.spikein-5.seq      >   $outdir/spikein-5.seq


outdir="FTO+/1-seq"
mkdir -p  $outdir
cat  Af_S15/Af_S15.spikein-1.seq      >   $outdir/spikein-1.seq
cat  Af_S15/Af_S15.spikein-2.seq      >   $outdir/spikein-2.seq
cat  Af_S15/Af_S15.spikein-3.seq      >   $outdir/spikein-3.seq
cat  Af_S15/Af_S15.spikein-4.seq      >   $outdir/spikein-4.seq
cat  Af_S15/Af_S15.spikein-5.seq      >   $outdir/spikein-5.seq



