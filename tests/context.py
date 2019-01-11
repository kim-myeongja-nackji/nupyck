import sys, os

pkg_path = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__), 
        '..'
    )
)

sys.path.append(pkg_path)
