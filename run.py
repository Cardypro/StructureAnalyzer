from pymol import cmd
import StructureAnalyzer_clean as sa

codes = ["1DVX","1NR6","1PXX","1SV9","2B17","2WEK","3CFQ","3IB0","3N8Y","4OJ4","4UBS","4XTA","4Z69","4ZBQ","4ZBR","5DBY","5U1R","6HN0","6HN1"]

for i in codes:
	sa.StructureAnalyzer(i,"DIF", autocompile = True)
	cmd.reinitialize()