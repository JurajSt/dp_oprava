# -*- coding: utf-8 -*-
# pripojenie k db
hostname = 'localhost'
username = 'postgres'
password = 'sql'
database = 'dp'
DB = "host=" + hostname + " user=" + username + " password=" + password + " dbname=" + database
# cesty k suborom
#nav_path = 'D:/diplomka/dp_projekt/data/test/navig/'#'D:/diplomka/dp_projekt/data/bordel/test/navig/'
obs_path = 'D:/diplomka/dp_projekt/data/kame_2/'  # 'D:/diplomka/telg_2/'
eph_path = 'D:/diplomka/dp_projekt/data/eph/eph/'  # igr18644.sp3'
csvpath = 'D:/diplomka/dp_projekt/data/SHMU/KAME_1617.csv'

#datum zaciatku a konca vypoctu
rok_start, mes_start, den_start = 2017, 1, 1
rok_konec, mes_konec, den_konec = 2017, 1, 22

satellite = ['G28']

# cisla druzic pre ktore bol v testovacich observacnych datach prijimany signal na L5, vychadza to z http://gpsworld.com/the-almanac/druzice block IIF
satellites_S5 = ['G01', 'G03','G06', 'G08', 'G09', 'G10', 'G24','G25', 'G26','G27', 'G30', 'G32']
# udaje o referencnej stanici
#nazovStanice = [[x,y,z], [elev. uhol od-do], [azimut od-do], referencna_vyska_stanice, nazov_stanice_rinex, nazov_meteostanice]
# pozn. suradnice vo WGS84, azimuta elev. uhol v stupnoch. Ak neznami azimut tak prazdny seznam,
# azimut zadavany podla wgs84 web mercator!!!
suradnice_vstup = 'xyz' # alebo 'blh'  tj zemepisne suradnice
kame = [[3892532.358, 1572220.333, 4785952.565], [5,25], [50, 80], 1.6, 'kame', 'Kamenica nad Cirochou']
telg = [[3947396.223, 1451396.020, 4780197.834], [5,30], [190, 220], 1.6, 'telg', 'Telgart']
hofn = [[2679690.298, -727951.336, 5722789.244], [5,25], [110, 150], 4.45, 'hofn', 'Hoefn']
ganp = [[3929181.9, 1455236.5, 4793653.8]]

pocet_observacii = 70  # minimalny pocet observacii pre analyzu

stanice = [kame]