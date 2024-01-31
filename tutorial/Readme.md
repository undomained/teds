Run these from the root source directory.

cd <teds>
export PYTHONPATH=$PWD/..:$PYTHONPATH
python GM/gm.py ../tutorial/gm.yaml
python SGM/sgm.py ../tutorial/sgm.yaml
./IM/tango_ckd_model/build/ckdmodel ../tutorial/im.cfg
./L1AL1B/tango_l1b/build/tango_l1b ../tutorial/l1al1b.cfg
python L1L2/l1bl2.py ../tutorial/l1l2.yaml