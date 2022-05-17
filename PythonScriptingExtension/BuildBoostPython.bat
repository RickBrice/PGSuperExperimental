REM copy ".\Supporting Documents\user-config-debugging-python.jam" %BOOSTDIR%
echo using python : %PYTHONVERSION% : %PYTHONDIR%\\python_d.exe : %PYTHONDIR%\\include : %PYTHONDIR%\\libs; > %BOOSTDIR%\user-config-debugging-python.jam
cd %BOOSTDIR%
call .\bootstrap.bat
call .\bjam --with-python python-debugging=off threading=multi variant=release link=shared address-model=64 stage
call .\bjam --user-config=user-config-debugging-python.jam --with-python python-debugging=on threading=multi variant=debug link=shared address-model=64 stage