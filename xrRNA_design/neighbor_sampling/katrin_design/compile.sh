clear &&
cd ../../ &&
python3 setup.py build && cp build/lib.linux-x86_64-cpython-312/infrared/ Test/ -r && cp build/lib.linux-x86_64-cpython-312/infrared/ Test/katrin_design -r;
cd Test/katrin_design
