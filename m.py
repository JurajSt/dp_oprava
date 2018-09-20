# -*- coding: utf-8 -*-
import datetime, time
from numpy import *
from input import *
import model
from azi_ele_coor import *
import matplotlib.pyplot as plt
from dataSHMU import *
import reader
from astropy.stats import LombScargle
#import poloha_druzice

ephfiles = np.sort(next(os.walk(eph_path))[2])

# datum zaciatku a konca vypoctu
dat_start = date(rok_start, mes_start, den_start)
dat_end = date(rok_konec, mes_konec, den_konec)
delta = dat_end-dat_start
pocet_dni = int((delta.total_seconds())/86400)

row = 1
col = 0
c = 0

doy_list1 = []
doy_list2 = []
SHMUdata = [[], []]
vyska1 = []
vyska2 = []
vyska3 = []

while c < len(stanice):

    SHMU = fSHMU(c) # spracovanie snehovych dat

    for f in ephfiles:

        if 'sp3' != f.split('.')[1]:
            continue

        e = eph_path + f # cesta k SP3 suboru


        sp3, date_e = reader.fReadSP3(e)  # nacitanie presnych efemerid
        doy_eph = reader.fdayofyear(int(date_e[0]), int(date_e[1]), int(date_e[2]))  # prepocet datumu na den v roku
        dat = date(int(date_e[0]), int(date_e[1]), int(date_e[2]))  # datum prveho zaznamu sp3 suboru
        # print int(date_e[0]), int(date_e[1]), int(date_e[2])

        if dat > dat_end:  # kontrola ci sa datum zaznamu sa nachadza v zadanom intervale datumu
            break
        if dat_start > dat:
            continue

        if 0 < doy_eph < 10:  # zistenie nazvu pre observacne data pomocou doy
            namefile = '00' + str(doy_eph) + '0'
        elif 10 <= doy_eph < 100:
            namefile = '0' + str(doy_eph) + '0'
        else:
            namefile = str(doy_eph) + '0'

        print
        print e

        rok = date_e[0][2:]
        obsfile = obs_path + stanice[c][4] + namefile + '.' + rok + 'o'
        if not os.path.exists(obsfile):
            print 'ziadny observacny subor'
            continue
        # else: continue
        r_obs = open(obsfile, "r")
        lines = r_obs.readlines()
        r_obs.close()

        print obsfile
        print dat

        # parsrovanie observacnych dat
        header, version, headlines, headlength, obstimes, sats, numallsvs, numsvs, date_o = reader.scan(lines)
        data, obstypes = reader.processBlocks(lines, header, obstimes, headlines, headlength, sats, numallsvs, satellite)

        del lines

        indexS1 = obstypes.index('L1') + 2
        indexS2 = obstypes.index('L2') + 2


        den_hodiny = []
        for sat in satellite:
            print ("pocitam pre: " + sat)
            t = []
            signals1 = []
            signals2 = []
            signals5 = []
            sinele = []
            a = []

            for k, d in enumerate(data):
                # print d
                # print k, sats[k]

                if sat not in sats[k]:
                    continue
                # print "pocitam pre: " + sat, k

                for i in d:  # jedna sekunda pre x druzic
                    # print i
                    sv = i[0]
                    if sat not in sv:
                        continue

                    lg = reader.fInterpolateSatelliteXYZ(sp3, sv, i[1])  # vypocet lagrangeovej interpolacie
                    azi, ele = fComputeAziEle(stanice[c][0], [lg[0], lg[1], lg[2]])  # vypocet elevacneho uhla, azimutu

                    # print azi, ele
                    ele_min = 5
                    ele_max = 30
                    azi_min = 60
                    azi_max = 80


                    #print i[indexS1], i[indexS2]
                    #print

                    s1 = i[indexS1].split("_")
                    if len(s1) == 1:
                        s1 = i[indexS1][-2:]
                        s1 = float(s1)*6
                    else:
                        s1 = float(s1[1]) * 6
                    s2 = i[indexS2].split("_")
                    if len(s2) == 1:
                        s2 = i[indexS2][-2:]
                        s2 = float(s2)*6
                    else:
                        s2 = float(s2[1]) * 6

                    #print s1, s2

                    # filtrovanie vhodnych signalov podla azimutu a elevacneho uhla
                    if ele_min <= ele <= ele_max and azi_min <= azi <= azi_max:
                        t.append(i[1])
                        sinele.append(m.sin(m.radians(ele)))
                        signals1.append(m.pow(10, s1 / 20))
                        signals2.append(m.pow(10, s2 / 20))
                        # print ele
                        break

            if len(sinele) < pocet_observacii:  # filter poctu observacii
                print "pre ", sat, " malo observacii - vypocet neprebehne", len(sinele)
                continue

            # print max(sinele), len(signals1), len(signals2)

            sinelexx = [[], []]
            s1 = [[], []]
            s2 = [[], []]
            t2 = [[], []]
            sinele3 = [[], []]

            for l in range(len(sinele) - 1):
                if sinele[l] - sinele[l + 1] > 0:
                    sinelexx[0].append(sinele[l])

                elif sinele[l] - sinele[l + 1] < 0:
                    sinelexx[1].append(sinele[l])
            ind = 0
            for q in range(len(sinelexx)):
                if len(sinelexx[q]) != 1:
                    continue
                ind = sinele.index(sinelexx[q][0])

            if ind == 0:
                s1[0] = signals1
                s2[0] = signals2
                sinele3[0] = sinele
                t2[0] = t
            else:
                s1[0] = signals1[:ind + 1]
                s2[0] = signals2[:ind + 1]
                sinele3[0] = sinele[:ind + 1]
                t2[0] = t[:ind + 1]

                s1[1] = signals1[ind + 1:]
                s2[1] = signals2[ind + 1:]
                sinele3[1] = sinele[ind + 1:]
                t2[1] = t[ind + 1:]

            for w in range(len(sinelexx)):

                if len(sinele3[w]) < pocet_observacii:
                    print "pre ", sat, " malo observacii (" + str(len(sinele3[w])) + ") - vypocet neprebehne"
                    continue

                print "pre " + sat + " pocet observacii: " + str(len(sinele3[w]))

                SNR1 = np.array(s1[w])  # zoznam na np.array
                SNR2 = np.array(s2[w])
                x = np.array(sinele3[w])
                timenum = np.array(t2[w])

                alt_prob = []
                i = 0
                while i < 4:
                    alt_prob.append(i)
                    i = i + 0.1
                # print alt_prob

                f = []
                for i in range(len(alt_prob)):  # prepocet vysky na frekvencu pre dlzky vln L1, L2, L5
                    f.append((2 * alt_prob[i]) / 0.190)
                freqS1 = np.array(f)

                f = []
                for i in range(len(alt_prob)):
                    f.append((2 * alt_prob[i]) / 0.244)
                freqS2 = np.array(f)

                deg = 6
                modelSNR1, residuals1, SNR1_vector = model.modelSNR(timenum, SNR1, deg)
                modelSNR1_1, residuals1_1, SNR1_vector_1 = model.modelSNR(x, SNR1_vector, deg)

                modelSNR2, residuals2, SNR2_vector = model.modelSNR(timenum, SNR2, deg)
                modelSNR2_1, residuals2_1, SNR2_vector_1 = model.modelSNR(x, SNR2_vector, deg)

                powerS1 = LombScargle(SNR1_vector_1, x).power(freqS1)
                powerS2 = LombScargle(SNR2_vector_1, x).power(freqS2)
                ind2S1 = np.ndarray.argmax(powerS1)
                ind2S2 = np.ndarray.argmax(powerS2)

                print alt_prob[ind2S1], alt_prob[ind2S2]  # , alt_prob[ind2S5]

                vyska1.append(alt_prob[ind2S1])
                doy_list1.append(dat)
                vyska2.append(alt_prob[ind2S2])

                if dat in SHMU[0] and dat in doy_list1:
                    sneh_ind = SHMU[0].index(dat)
                    SHMUdata[0].append(SHMU[0][sneh_ind])
                    sneh = (SHMU[1][sneh_ind] / 100.0)
                    SHMUdata[1].append(SHMU[1][sneh_ind])

    if len(doy_list1) == 0:
        print 'Ziadny vypocet vysky. Skontroluj data'
        os.system("pause")
        sys.exit()
    # priemer vysok
    mean1 = np.mean(vyska1)
    mean2 = np.mean(vyska2)

    # graf 1 vypocitane vysky a referencna vyska
    m1 = [stanice[c][3], stanice[c][3]]
    dm = [min(doy_list1), max(doy_list1)]
    plt.plot(doy_list1, vyska1, 'o', markersize=10, color="blue")
    plt.plot(doy_list1, vyska2, '^', markersize=10, color="red")

    plt.plot(dm, m1, color="black")
    #plt.savefig(cesta_vystup + stanice[c][4] + '_ref_vyska.png')
    plt.show()

    vyska_r1 = []
    for v in vyska1:
        vyska_r1.append((kame[3] - v) * 100)

    vyska_r2 = []
    for v in vyska2:
        vyska_r2.append((kame[3] - v) * 100)

        # vypocet RMSE pre jednotlive signaly
        print
        try:
            rmse_val = reader.rmse(np.array(vyska_r1), np.array(SHMUdata[1]))
            print("rms chyba L1 is: " + str(rmse_val))
            rmse_val = reader.rmse(np.array(vyska_r2), np.array(SHMUdata[1]))
            print("rms chyba L2 is: " + str(rmse_val))
        except ValueError as e:
            print "ValueError: " + e
            print "datum observacie nie je v datach v meranej vyske snehu"

        # graf pre odhadovanu vysku snehu
        plt.plot(doy_list1, vyska_r1, 'o', markersize=10, color="blue")
        plt.plot(doy_list1, vyska_r2, '^', markersize=10, color="red")
        plt.plot(SHMUdata[0], SHMUdata[1])
        plt.show()

    c += 1