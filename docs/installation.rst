Installation
============

To install NPlib, you have two options: cloning the repository or using pip.

Option 1: Clone the repository
------------------------------

1. **Clone the repository**:

    Open a terminal and run:
    
    ```
    git clone https://github.com/farrisric/NPlib.git
    ```

2. **Navigate to the project directory**:

    ```
    cd NPlib
    ```

3. **Create a virtual environment** (optional but recommended):

    ```
    python3 -m venv venv
    source venv/bin/activate
    ```

4. **Install the required dependencies**:

    ```
    pip install -r requirements.txt
    ```

5. **Install NPlib**:

    ```
    python setup.py install
    ```

6. **Verify the installation**:

    ```
    python -c "import nplib; print(nplib.__version__)"
    ```

Option 2: Install via pip
-------------------------

1. **Install NPlib directly from PyPI**:

    ```
    pip install npl
    ```

2. **Verify the installation**:

    ```
    python -c "import nplib; print(nplib.__version__)"
    ```

If you encounter any issues during installation, please refer to the `README.md` file or open an issue on the [GitHub repository](https://github.com/farrisric/NPlib/issues).