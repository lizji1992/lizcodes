hmr=~/panfs/bin/methpipe/bin/hmr

methdir=$1
outdir=$2

for f in $(find -L $methdir -name "*.meth")
do
    echo runing hmr on $f
    fname=$(basename $f)
    outf=${fname/.meth/.hmr}
    [ ! -f $outdir/$outf ] \
        && $hmr -o $outdir/$outf $methdir/$fname
    echo hmr finished on $f
done

