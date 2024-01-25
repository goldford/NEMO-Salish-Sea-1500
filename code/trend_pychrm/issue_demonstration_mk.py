# I don't understand why I get different trends per season
import numpy as np
import mannkendall as mk
import datetime as dt

print(mk.__version__)
import sys; print('Python %s on %s' % (sys.version, sys.platform))
# 36 years of fictional seasonal anomalies
winter_dat = np.array([-0.23001771, -0.84052181, -1.23213363, -0.76677547, -0.37802406,
       -1.00168761, -0.19810082, -0.00479787, -0.3879219 ,  0.1598607 ,
        0.21215476,  0.61506259, -0.24934041,  0.60148242,  0.26912125,
       -0.41568088,  0.16919191,  0.33381554,  0.19711335, -0.3606954 ,
        0.10904288,  0.01487721,  0.36594954,  0.05137299,  0.52846712,
        0.37628641,  0.14614995, -0.01926223,  0.56417156,  0.27049467,
        0.20427819,  0.19442583,  0.35200845,         np.nan,  0.08953783,
        0.60767666])

spring_dat = np.array([-0.06035071, -0.988043  , -1.39324562, -0.50995384, -0.85708491,
       -0.80519989, -0.3801143 , -0.08300517, -0.10336021, -0.37871942,
       -0.13816126,  0.44140764, -0.25222606,  0.46402375,  0.37023428,
        0.03578017, -0.41652947,  0.27963351,  0.20860907, -0.00739693,
        0.076229  , -0.06976548,  0.80467958,  0.31614347,  0.53311702,
        0.26518428,  0.27314458, -0.19927579,  0.71170745,  0.00906855,
        0.15378613,  0.14391438, -0.27978446,         np.nan,  0.88902366,
        0.90677563])

summer_dat = np.array([ 0.22557246, -0.71489365, -1.32649691, -0.90483987, -0.92756592,
       -0.76139591, -0.33367437, -0.50860319, -1.04543393, -0.47956156,
       -0.47048985, -0.28036398, -0.08423861,  0.41584444,  0.12030349,
        0.28267624, -0.1489518 ,  0.18771227,  0.04758116,  0.09951411,
        0.38648407, -0.02137743,  0.6961336 ,  0.22047614,  0.5088918 ,
        0.28604842,  0.12163195, -0.02130959,  1.06345857,  0.1090539 ,
        0.41110131,  0.36717935,  0.3645205 ,         np.nan,  0.89566435,
        1.2509044 ])
fall_dat = np.array([-0.49176537, -0.89438428, -1.22146029, -0.66709172, -0.36386329,
       -0.744658  , -0.44863015, -0.3654154 ,  0.18158329, -0.06274773,
        0.05225475,  0.07571491, -0.0948965 ,  0.34798375, -0.16340824,
       -0.29530515,  0.04672506,  0.33498664, -0.010352  , -0.15621559,
        0.22566208,  0.13370434,  0.48756215,  0.4523201 ,  0.41637158,
        0.31148928,  0.25684139,  0.49534467,  0.52922754,  0.19323551,
        0.14456546,  0.19001399, -0.01786903,         np.nan,  0.7866529 ,
        0.38562555])

# REDO DATES
# original code in demo example has error where first two date elements are same year
single_seas_dts = np.array([dt.datetime(yr, 1, 1) for yr in range(1970,2006,1)])
print(single_seas_dts)

winter_seas_dts = single_seas_dts
spring_seas_dts = np.array([dt.datetime(yr, 4, 1) for yr in range(1970,2006,1)])
summer_seas_dts = np.array([dt.datetime(yr, 7, 1) for yr in range(1970,2006,1)])
fall_seas_dts = np.array([dt.datetime(yr, 10, 1) for yr in range(1970,2006,1)])

print("Test to verify indiv seas trends match when combined together")
resolution = 0.001
print("== Seasons Analysed Individually ==")
print("winter")
print(mk.mk_temp_aggr(winter_seas_dts,winter_dat, resolution))
print("spring")
print(mk.mk_temp_aggr(spring_seas_dts,spring_dat, resolution))
print("summer")
print(mk.mk_temp_aggr(summer_seas_dts,summer_dat, resolution))
print("fall")
print(mk.mk_temp_aggr(fall_seas_dts,fall_dat, resolution))

print("== Seasons Analysed In One Dataset ==")
multi_seas_dts = [winter_seas_dts,spring_seas_dts, summer_seas_dts, fall_seas_dts]
multi_seas_data = [winter_dat, spring_dat, summer_dat, fall_dat]
out = mk.mk_temp_aggr(multi_seas_dts,multi_seas_data, resolution)

for n in range(len(out)-1):
    print("Season {}: {}".format(n+1,out[n]))
print("Overall trends: ", out[len(out)-1])

# rand_mult = 0.01
# fake_data = 0.5 * np.ones((36,)) + (np.arange(36)*0.01) + (np.random.rand(36)*rand_mult)
# fake_data = 0.5 * np.ones((36,)) + (np.arange(36)*0.01)

# print("all seasons included, same dates:")
# all_seasons_dat = [winter_dat, spring_dat, summer_dat, winter_dat]
# dates_seasons = [dates_ar, dates_ar, dates_ar, dates_ar]
# out = mk.mk_temp_aggr(dates_seasons,all_seasons_dat, 0.0001)
# print(out)
