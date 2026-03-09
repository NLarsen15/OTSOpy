#!/bin/bash
# Simple script to build MiddleMan extension with f2py
# This is the working approach that generates proper PyInit_MiddleMan symbols

echo "Building MiddleMan extension with f2py..."
python -m numpy.f2py -c --fcompiler=gnu95 -m MiddleMan *.f *.f95 *.for

if [ $? -eq 0 ]; then
    echo "Build successful!"
    
    # Find the generated extension file
    EXTFILE=$(ls MiddleMan.cpython-*.so 2>/dev/null || ls MiddleMan.so 2>/dev/null)
    if [ -n "$EXTFILE" ]; then
        echo "Generated extension: $EXTFILE"
        
        # Verify PyInit_MiddleMan symbol exists
        if nm -D "$EXTFILE" | grep -q PyInit_MiddleMan; then
            echo "✓ Extension has proper PyInit_MiddleMan symbol"
            
            # Test import
            echo "Testing import..."
            python -c "import MiddleMan; print('✓ MiddleMan imports successfully')" 2>/dev/null
            if [ $? -eq 0 ]; then
                echo "✓ Extension works correctly!"
            else
                echo "✗ Extension import failed"
            fi
        else
            echo "✗ Extension missing PyInit_MiddleMan symbol"
        fi
    else
        echo "✗ No extension file found"
    fi
else
    echo "Build failed!"
    exit 1
fi