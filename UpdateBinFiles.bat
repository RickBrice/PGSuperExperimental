REM - Script to prepare for Release

SET BINTARGET=bin
SET REGFREECOM=\ARP\BridgeLink\RegFreeCOM


REM - Experimental Extensions
xcopy /Y /d %REGFREECOM%\x64\Release\IfcExtensions.dll	%BINTARGET%\Extensions\Experimental\x64\
