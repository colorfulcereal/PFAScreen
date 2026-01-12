#!/bin/bash
# Installation script for PFAScreen on macOS
echo "-----------------------------------------"
echo ""
echo "PFAScreen Installation for macOS"
echo ""
echo "Checking Python installation..."

# Check if Python 3.9 is installed
if command -v python3 &> /dev/null; then
    PYTHON_VERSION=$(python3 --version 2>&1 | awk '{print $2}')
    echo "Python version found: $PYTHON_VERSION"

    # Check if version is 3.9 or higher
    MAJOR=$(echo $PYTHON_VERSION | cut -d. -f1)
    MINOR=$(echo $PYTHON_VERSION | cut -d. -f2)

    if [ "$MAJOR" -eq 3 ] && [ "$MINOR" -ge 9 ]; then
        echo "Python version is compatible!"
    else
        echo "Warning: Python 3.9 or higher is recommended."
        echo "Current version: $PYTHON_VERSION"
        echo "Please install Python 3.9+ from https://www.python.org/downloads/"
        read -p "Press Enter to continue anyway or Ctrl+C to exit..."
    fi
else
    echo "Python 3 not found!"
    echo "Please install Python 3.9+ from https://www.python.org/downloads/"
    read -p "Press Enter to continue or Ctrl+C to exit..."
    exit 1
fi

echo ""
echo "Installing required Python packages..."
echo ""

# Ensure pip is installed
python3 -m ensurepip --default-pip 2>/dev/null || true

# Upgrade pip
python3 -m pip install --upgrade pip

# Install requirements
python3 -m pip install -r requirements.txt

if [ $? -eq 0 ]; then
    echo ""
    echo "Installation successfully finished!"
    echo "You can now run PFAScreen using: ./run_pfascreen.sh"
    echo ""
else
    echo ""
    echo "Installation encountered errors. Please check the output above."
    echo ""
fi

read -p "Press Enter to exit..."
