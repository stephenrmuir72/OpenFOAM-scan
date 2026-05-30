# ----------------------------------------------------------------------------
# Common Mingw cross-compile preferences
# ----------------------------------------------------------------------------

export FOAM_CONFIG_ETC="etc-mingw"
export WM_COMPILER=Mingw

# Detect host architecture
ARCH=$(uname -m)

if [ "$ARCH" = "aarch64" ]; then
    # ---------------- ARM64 HOST ----------------
    export WM_ARCH_OPTION=ARM64
    export WM_CROSS_COMPILE=aarch64-w64-mingw32
    export WM_OPTIONS=winARM64Mingw
    export WM_MPLIB=MSMPI

elif [ "$ARCH" = "x86_64" ]; then
    # ---------------- x86_64 HOST ----------------
    unset WM_ARCH_OPTION
    unset WM_CROSS_COMPILE
    export WM_OPTIONS=win64Mingw
    export WM_MPLIB=msmpi-10.0

else
    echo "Unsupported architecture: $ARCH"
    return 1 2>/dev/null || exit 1
fi

# Common cleanup
unset WM_COMPILE_CONTROL