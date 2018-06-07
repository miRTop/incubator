set -o pipefail  # trace ERR through pipes
set -o errtrace  # trace ERR through 'time command' and other functions
set -o nounset   ## set -u : exit the script if you try to use an uninitialised variable
set -o errexit   ## set -e : exit the script if any statement returns a non-true return value
set -v
set -x
export PS4='+(${BASH_SOURCE}:${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

ROOT=`pwd`

cd tewari

cd bcbio
mirtop gff --hairpin ../../reference/hairpin.fa --gtf ../../reference/mirbase.gff3  -o mirtop --sps hsa --format seqbuster */*mirbase-ready.counts
mirtop stats -o stats mirtop/*ready.gtf
mirtop counts --gff mirtop/mirtop.gtf --out mirtop --hairpin ../../reference/hairpin.fa --gtf ../../reference/mirbase.gff3 --add-extra --sps hsa
cd ..

cd mirge
mirtop stats -o stats *gff 
mirtop gff --hairpin ../../reference/hairpin.fa --gtf ../../reference/mirbase.gff3  -o mirtop --sps hsa --format gff *gff
mirtop counts --gff mirtop/mirtop.gtf --out mirtop --hairpin ../../reference/hairpin.fa --gtf ../../reference/mirbase.gff3 --add-extra --sps hsa
# samples have a different names than files, get the relationship
# to normalize steps across tools
python $ROOT/scripts/get_sample_names_from_gff.py *gff > sample_fn_key.txt
cd ..

cd isomirsea
	mirtop  gff --hairpin ../../reference/hairpin.fa --gtf ../../reference/mirbase.gff3  -o mirtop --sps hsa --format isomirsea mirbase21/*gff
	mirtop stats -o mirtop mirtop/*all.gtf
	mirtop counts --gff mirtop/mirtop.gtf --out mirtop --hairpin ../../reference/hairpin.fa --gtf ../../reference/mirbase.gff3 --add-extra --sps hsa
cd ..

cd srnabench
    for folder in `ls .`; do echo $folder;  mirtop  gff --hairpin ../../reference/hairpin.fa --gtf ../../reference/mirbase.gff3  -o $folder --sps hsa --format srnabench $folder; done
    mirtop stats -o mirtop */mirtop.gtf
    mirtop gff --hairpin ../../reference/hairpin.fa --gtf ../../reference/mirbase.gff3  -o mirtop --sps hsa --format gff */mirtop.gtf
    mirtop counts --gff mirtop/mirtop.gtf --out mirtop --hairpin ../../reference/hairpin.fa --gtf ../../reference/mirbase.gff3 --add-extra --sps hsa
cd ..
