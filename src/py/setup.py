from setuptools import setup

setup(
    name="gamma",
    version="1.0",
    description="Gamma pyf90 interface package",
    url="/",
    author="INRAE/Hydris Hydrologie",
    packages=["gamma","gamma.core","gamma.wrapping", "gamma.smashplug"],
    install_requires=[
    ],
    package_data={
        "gamma": ["gamma/wrapping/_wrapping.cpython-38-x86_64-linux-gnu.so","gamma/wrapping/_wrapping.cpython-310-x86_64-linux-gnu.so"]
    },
    include_package_data=True,
    zip_safe=False,
)
