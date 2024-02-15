source /nobackup/users/leune/venv/tango/bin/activate

cd ../../teds
export PYTHONPATH=$PWD/..:$PYTHONPATH

python ./GM/gm.py ../cfg/nitro/gm.yaml
python ./SGM/sgm_no2.py ../cfg/nitro/sgm.yaml