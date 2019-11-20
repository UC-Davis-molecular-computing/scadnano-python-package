source ../../scadnano-dart/define_root_dir.sh

DOC_SUBDIR="docs/"
DOC_DIR="$ROOT_DIR/$DOC_SUBDIR"

echo Assuming root directory of website is $ROOT_DIR
echo "uploading index.html first to $DOC_DIR"

# upload this first since it is the main document and is fast
scp _build/html/index.html $DOC_DIR

echo "now uploading the rest to $DOC_DIR"

# now upload everything else
scp -r _build/html/* $DOC_DIR