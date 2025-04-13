(Get-Content run_fit.py) | foreach { $_ -replace 'sname = ""','sname = "PGM"' } | Set-Content run_fit.py
python run_fit.py
(Get-Content run_fit.py) | foreach { $_ -replace 'sname = "PGM"','sname = "GaPGM"' } | Set-Content run_fit.py
python run_fit.py
(Get-Content run_fit.py) | foreach { $_ -replace 'sname = "GaPGM"','sname = "InPGM"' } | Set-Content run_fit.py
python run_fit.py
(Get-Content run_fit.py) | foreach { $_ -replace 'sname = "InPGM"','sname = "SnPGM"' } | Set-Content run_fit.py
python run_fit.py
(Get-Content run_fit.py) | foreach { $_ -replace 'sname = "SnPGM"','sname = ""' } | Set-Content run_fit.py