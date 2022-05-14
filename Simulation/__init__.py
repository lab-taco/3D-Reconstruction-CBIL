""" OS Check """
import os
import platform
import sysconfig
import sys
print('OS Check...')
print("os.name                     ", os.name)
print("sys.platform                ", sys.platform)
print("platform.system()           ", platform.system())
print("sysconfig.get_platform()    ", sysconfig.get_platform())
print("platform.machine()          ", platform.machine())
print("platform.architecture()     ", platform.architecture())
if platform.system()!= "Windows":
    import txaio
    txaio.use_asyncio()
""" -------------- """