source /nobackup/users/leune/venv/tango/bin/activate

cd ../teds
export PYTHONPATH=$PWD/..:$PYTHONPATH

python ./GM/gm.py ../no2/gm_no2.yaml
python ./SGM/sgm_no2.py ../no2/sgm_no2.yaml