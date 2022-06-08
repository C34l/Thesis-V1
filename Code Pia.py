"""Minimum der FSQ der Temperaturen wird gesucht. Neben der experimentellen Temperatur wird die berechnete Temperatur
 aus der SLE-Gleichung bestimmt. Zur Minimumberechnung schon in Python verfügbare Fkt. genutzt.

Werte aus jr9430000231.pdf

xA0 = xS und xB0 =xWasser

! bei Berechnung der Temp zur Modellierung noch die 100 Datenpunkte(Zusammensetzungen) entsprechend x_experiment anpassen"""
import pandas as pd
import scipy.optimize as opt
import math
import numpy as np

# Konstanten:
T0S_LS = 404.75  # Kelvin, bzw. 131,65 °C
deltahLS_R = 24500.0  # J/mol
R = 8.31446261815324  # kg * m ^ 2 / (s ^ 2 * mol * K)

# Variabeln-Listen
T0 = [297.65, 300.65, 304.65, 310.15, 314.65, 317.15, 318.65, 321.65, 323.65, 325.65, 327.65, 330.15, 333.65, 341.15]  # Kelvin
xA0 = [0.010425105, 0.012737984, 0.016645163, 0.024029341, 0.034293222, 0.044991658, 0.061853946, 0.080539095, 0.10051017, 0.117793448, 0.132558662, 0.155114202, 0.18675243, 0.255024557]
# T0 = [389.438, 391.739, 393.963, 394.04, 394.193, 399.101, 401.632, 404.546]
# xA0 = [0.702997, 0.747956, 0.798365, 0.811989, 0.813352, 0.899183, 0.941417, 1.00272]
xB0 = []
for i in range(len(xA0)):
    xB0.append(1-xA0[i])

# Startwerte für Iteration
delta_g_0 = np.array([-2842, 9979])

#dictionary for all values
arg_temp = {'deltahLS_R': deltahLS_R,
            'T0S_LS': T0S_LS,
            'xA0': xA0,
            'xB0': xB0,
            'T0': T0,
            'R': R,
            'delta_g_0': delta_g_0}     # so schreiben, wenn im Nachhinein noch etwas hinzufügen: arg_temp['delta_g_0'] = delta_g_0


def calc_gamma_mix(x0, i, Ti_mod, **arg_temp):
    R = arg_temp['R']
    xA0 = arg_temp['xA0']
    xB0 = arg_temp['xB0']

    g_AB = x0[0]
    g_BA = x0[1]

    global alpha
    alpha = 0.4 # alphaAB =alphaBA
    Ti_mod = Ti_mod
    tau_BA = g_BA/(R*Ti_mod)
    G_BA = math.exp(-alpha*tau_BA)
    tau_AB = g_AB/(R*Ti_mod)
    G_AB = math.exp(-alpha*tau_AB)

    # print(tau_AB)
    # print(tau_BA)
    # print(G_AB)
    # print(G_BA)

    gamma = math.exp(xB0[i] ** 2 * ((tau_BA * (G_BA / (xA0[i] + xB0[i] * G_BA)) ** 2 + tau_AB * G_AB / (xA0[i] * G_AB + xB0[i]) ** 2)))

    return gamma


def Temp_calc(T, x0, i, deltahLS_R, T0S_LS, R):   # funktioniert mit arg_temo irgendwie nicht **arg_temp):  #maybe better to calc only one point instead of xA0, for also using this def later while plotting, without need to overwrite xA0

    Ti_mod = T[0]  # wie sonst definieren? Eigentlich ja in x0, aber da doch delta_gij (Parameter anzupassen) drin?
    gamma = calc_gamma_mix(x0, i, Ti_mod, **arg_temp)  # da jedes Mal neu calc_gamma aufgerufen, array immer neu überschrieben
    f = 1 / (gamma * xA0[i]) * np.exp(-deltahLS_R / (R * Ti_mod) * (1 - Ti_mod / T0S_LS)) - 1   # umstellen der SLE-Gleichung nach 0 = ... dann fsolve anwenden
    #f = (-deltahLS_R / (R * Ti_mod) * (1 - Ti_mod / T0S_LS)) / np.log(gamma * xA0[i]) - 1
    #f = np.log(gamma * xA0[i])/(-deltahLS_R / (R * Ti_mod) * (1 - Ti_mod / T0S_LS)) - 1
    return f


def min_FQS(x0, **arg_ls):  # dictionary arg_temp, but other name in this function, gets dictionary by calling function in ls(least_squares)
    xA0 = arg_temp['xA0']
    T0 = arg_temp['T0']
    deltahLS_R = arg_temp['deltahLS_R']
    T0S_LS = arg_temp['T0S_LS']
    R = arg_temp['R']

    #Ti_mod_0 = 273.15   # Ti_mod_0 beschreibt Startwert bei der Suche
    global Ti_mod
    Ti_mod = np.zeros(14)
    for i in range(len(xA0)):
        Ti_mod_0 = T0[i]
        Ti_mod[i] = opt.fsolve(Temp_calc, Ti_mod_0, args=(x0, i, deltahLS_R, T0S_LS, R))     # fsolve(Funktion, Grenzen), "args" nutzen, um Extra-Argumente zu übergeben #Ti_mod_0 ist Variable(oben T), die verändert werden soll
        if abs(Ti_mod[i]-T0[i]) < 1e-6:
            print("Fehler!")

    global FQS_Terme  #damit auch außerhalb aufgerufen werden kann
    FQS_Terme = (T0[0:len(xA0)] - Ti_mod) ** 2  #* 10**6
    global normFQS_Terme
    normFQS_Terme = ((T0[0:len(xA0)] - Ti_mod) / T0[0:len(xA0)]) ** 2
    # FQS_Terme = ((T0[0:len(xA0)] - Ti_mod)/T0[0:len(xA0)])**2   # normiert

    global Obj
    Obj = np.sum(FQS_Terme)
    #LS = (T0[0:len(xA0)] - Ti_mod)
    return Obj


matrix_comparison = [] #np.zeros([80, 80])
# matrix_comparison.append(xA0)
# matrix_comparison.append(T0)

matrix_comparison.append(delta_g_0)
res_1 = opt.minimize(min_FQS, delta_g_0, method='Powell')  #Nelder-Mead')   # args={**arg_temp})

print('ARD [%] = ' + str(100 * np.sum(np.abs(np.divide(-res_1.fun, T0)))))   # res_1.fun calls function min_FQS at A1 & A2, better than calling the function again, because value is already calculated
print('deltag_AB [ ] = ' + str(res_1.x[0]))   # str dazu, wenn + (da float und string nicht mgl.)
print('deltag_BA [ ] = ' + str(res_1.x[1]))
print('alpha =', alpha)

x0 = [res_1.x[0], res_1.x[1]]
T_calc = Ti_mod     # wird trotzdem ohne Fehler geprinted... Was ist das Problem?

T_diff_norm = np.zeros(14)
T_diff = np.zeros(14)
for i in range(len(T0)):
    T_diff_norm[i] = abs((T0[i]-T_calc[i])/T0[i])
    T_diff[i] = abs((T0[i]-T_calc[i]))
ARD_neu_norm = 100/len(T0) * sum(T_diff_norm)
ARD_neu = 100/len(T0) * sum(T_diff)
print('ARD_neu_norm_bzgl.T [%] =', ARD_neu_norm)
print('ARD_neu_bzgl.T [%] =', ARD_neu)

FQS = FQS_Terme #min_FQS(x0, **arg_temp) # jetzt Summe davon als Rückgabe der min_function --> vorher Vektor global speichern
print("Die einzelnen Fehlerquadratsummen lauten:", FQS)
# print("Zu den Molenbrüchen xA0:", xA0)
# print("Ti_mod sind gleich", T_calc)
# T_calc = np.zeros(14)     #ab hier iwas falsch 4 Zeilen--> liefert teilweise richtige, teilweise falsche Temperaturen
# for i in range(len(xA0)):
#     T_calc[i] = T0[i] - (FQS[i])**0.5  #/10**6
    # T_calc[i] = (1-FQS[i]**0.5) * T0[i]   # eigentlich müsste das die Formel sein!! FQS zwar gleich klein, aber dann kommen falsche Temperaturen raus :/passt nicht zur Kurve des Experiments

print('\n', 'xA0     ', 'T_exp   ', 'T_calc')
for i in range(len(xA0)):
    print(round(xA0[i], 4), '   ', round(xB0[i], 4), '  ', T0[i], ' ', T_calc[i].round(15))  # would like to save all data in array - but how?, printing array directly would look much better


matrix_comparison.append(xA0)
matrix_comparison.append(T0)
matrix_comparison.append(T_calc)

print(np.sum(res_1.fun))
print("Die (im Moment) nicht normierte Fehlerquadratsumme beträgt:", min_FQS(x0, **arg_temp))   #hier jetzt Summe direkt in min_FQS, also ncihtmehr minimize_FQS aufrufen
print("Die nicht normierte Fehlerquadratsumme beträgt:", sum(FQS_Terme))
print("Die normierte Fehlerquadratsumme beträgt:", sum(normFQS_Terme))

# matrix_comparison = np.zeros((9, 3))
# matrix_comparison[:, 0] = xA0   # alle Stellen, 1. Spalte
# matrix_comparison[:, 1] = T0
# matrix_comparison[:, 2] = T_calc
# #matrix_comparison[:, 3] = gamma_neu

matrix = pd.DataFrame(matrix_comparison)
matrix.to_excel("Minimize-S(undH2O).xlsx", sheet_name='S(undH2O)-Löslichkeit-NRTL')

# np.set_printoptions(precision=5)  # ab hier werden alle nachfolgenden arrays auf 5 Nachkommatsellen gerundet
# print("\n", "xA0:        T_exp:       T_calc: ")
# print(matrix_comparison)

for i in range(len(xA0)):
    print("Der Aktivitätskoeffizient für xA0[", i, "]", calc_gamma_mix(x0, i, T_calc[i], **arg_temp))
# print("Oder diret so, die nicht normiete FSQ über globale Variable ist: ", Obj)



points = 100     # following gives 100 points, using this for plotting later
xstart = 0.0001
xfinish = 0.27  #xA0[-1]     # xA0 --> Löslichkeiten und -1 heißt letztes Element # 1- x0 stand hier, aber ergibt doch keinen Sinn oder?
xA0_neu = np.linspace(xstart, xfinish, points)
arg_temp['xA0_neu'] = xA0_neu  #overwriting xA0, 1000 instead of 9 points for plotting ! write AFTER using xA0

xB0_neu = np.zeros(100)
for i in range(0, len(xA0_neu)):
    xB0_neu[i] = (1-xA0_neu[i])
arg_temp['xB0_neu'] = xB0_neu


def curve_fit_Temp (Ti_mod, i, deltag, xA0, xB0, deltahLS_R, T0S_LS, R):    # **arg_ls):
    g_AB = deltag[0]
    g_BA = deltag[1]
    alpha = 0.4  # alphaAB =alphaBA
    tau_BA = g_BA / (R * Ti_mod)
    G_BA = math.exp(-alpha * tau_BA)
    tau_AB = g_AB / (R * Ti_mod)
    G_AB = math.exp(-alpha * tau_AB)

    # print(xA0[i])
    # funcneu = []
    # for Ti_mod in T_modell_0:
    #!!!! hier jetzt mit xB0 rechnen
    funcneu = (xB0[i] ** 2) * ((tau_BA * (G_BA / (xA0[i] + xB0[i] * G_BA)) ** 2 + tau_AB * G_AB / (xA0[i] * G_AB + xB0[i]) ** 2)) - math.log(1 / xA0[i] * math.exp(-deltahLS_R / (R * Ti_mod) * (1 - Ti_mod / T0S_LS)))
    return funcneu

T_modell = np.zeros(100)  #[]
T_modell_0 = np.zeros(100)  # K einfach Startwert für alle x beim Plotten

print('xA0    ', 'T_modell(Plot)')
for i in range(0, len(xA0_neu)):
    T_modell_0[i] = 300
    T_modell[i] = (opt.fsolve(curve_fit_Temp, T_modell_0[i], args=(i, x0, xA0_neu, xB0_neu, deltahLS_R, T0S_LS, R)))
    print(xA0_neu[i], '    ', T_modell[i])

# print("100 Temperaturen zur Modellierung der Kurve", T_modell)

WertexA0 = pd.DataFrame(xA0_neu)
WertexB0 = pd.DataFrame(xB0_neu)
WerteT_modell = pd.DataFrame(T_modell)

writer = pd.ExcelWriter('Python_Curve_fit_S(undH20).xlsx', enigine='xlsxwriter') #engine soll man einfach so angeben laut internet, Fehlermeldung ist nicht wichtig

WertexA0.to_excel(writer, sheet_name='xA0')
WertexB0.to_excel(writer, sheet_name='xB0')
WerteT_modell.to_excel(writer, sheet_name='T_modell')

writer.save()


# neuer Plott mit viel mehr Wertzen und viel größerem Bereich
points = 2000     # following gives 100 points, using this for plotting later
xstart = 0.0001
xfinish = 0.5  #xA0[-1]     # xA0 --> Löslichkeiten und -1 heißt letztes Element # 1- x0 stand hier, aber ergibt doch keinen Sinn oder?
xA0_neu = np.linspace(xstart, xfinish, points)
arg_temp['xA0_neu'] = xA0_neu  #overwriting xA0, 1000 instead of 9 points for plotting ! write AFTER using xA0

xB0_neu = np.zeros(2000)
for i in range(0, len(xA0_neu)):
    xB0_neu[i] = (1-xA0_neu[i])
arg_temp['xB0_neu'] = xB0_neu


def curve_fit_Temp (Ti_mod, i, deltag, xA0, xB0, deltahLS_R, T0S_LS, R):    # **arg_ls):
    g_AB = deltag[0]
    g_BA = deltag[1]
    alpha = 0.4  # alphaAB =alphaBA
    tau_BA = g_BA / (R * Ti_mod)
    G_BA = math.exp(-alpha * tau_BA)
    tau_AB = g_AB / (R * Ti_mod)
    G_AB = math.exp(-alpha * tau_AB)

    # print(xA0[i])
    # funcneu = []
    # for Ti_mod in T_modell_0:
    #!!!! hier jetzt mit xB0 rechnen
    funcneu = (xB0[i] ** 2) * ((tau_BA * (G_BA / (xA0[i] + xB0[i] * G_BA)) ** 2 + tau_AB * G_AB / (xA0[i] * G_AB + xB0[i]) ** 2)) - math.log(1 / xA0[i] * math.exp(-deltahLS_R / (R * Ti_mod) * (1 - Ti_mod / T0S_LS)))
    return funcneu

T_modell = np.zeros(2000)  #[]
T_modell_0 = np.zeros(2000)  # K einfach Startwert für alle x beim Plotten

print('xA0    ', 'T_modell(Plot)')
for i in range(0, len(xA0_neu)):
    T_modell_0[i] = 300
    T_modell[i] = (opt.fsolve(curve_fit_Temp, T_modell_0[i], args=(i, x0, xA0_neu, xB0_neu, deltahLS_R, T0S_LS, R)))
    print(xA0_neu[i], '    ', T_modell[i])

# print("100 Temperaturen zur Modellierung der Kurve", T_modell)

WertexA0 = pd.DataFrame(xA0_neu)
WertexB0 = pd.DataFrame(xB0_neu)
WerteT_modell = pd.DataFrame(T_modell)

writer = pd.ExcelWriter('Python_Curve_fit_S(undH20)_größererBereich.xlsx', enigine='xlsxwriter') #engine soll man einfach so angeben laut internet, Fehlermeldung ist nicht wichtig

WertexA0.to_excel(writer, sheet_name='xA0')
WertexB0.to_excel(writer, sheet_name='xB0')
WerteT_modell.to_excel(writer, sheet_name='T_modell')

writer.save()







# ab hier alles nur noch Extras

'''Jetzt nochmal so umschreiben/..., dass man manuell Parameter eingeben kann und die netsprechenden Temperaturen 
damit ausgerechnet werden (aber schon auch zu x_experiment) --> jetzt aber nicht die FSQ minimieren, sondern nur einmal 
Nullstellen für die Parameter für alle xA0 finden und dann das ausgeben --> keine Veränderung von deltagij/deltagji '''

deltag_AB = -2827.972479393138
deltag_BA = 9907.40489115058
deltag = [deltag_AB, deltag_BA]
print("Die manuell eingegebenen Parameter sind:")
print("deltag_AB =", deltag_AB)
print("deltag_BA", deltag_BA)

sum_FQS_manuell = min_FQS(deltag, **arg_temp)   #berechnet nur Fehlerquadratsumme, Minimiert diese noch nicht
print("Die Fehlerqudadratsumme der manuell eingegebenen Parametern *10^6 ist:")
print(sum_FQS_manuell)
print("Die Temp zu manuell eingegebenen Parametern sind", Ti_mod)   #ersetzt alle weiteren Zeilen, wenn einfach oben das array global gespeichert wird --> verlangsamt nur wahrscheinlich das Programm ziemlich

# # FQS_man =[]
# FQS_man = FQS_Terme   # geht, da globale Variable, aber !!!! FQS_Terme ist damit also überschrieben (da DokumentenENDE)
# T_calc_man = np.zeros(14)  # man = manuell
# for i in range(len(xA0)):
#     T_calc_man[i] = T0[i] - (FQS_man[i]) ** 0.5   #!!! /10**6, da vorher mal eine Million (weil Annäherung komisch)
#
# print("Die kalkulierten Temperaturen für die manuell eingegebenen Paramater lauten:")
# print(T_calc_man)



''' Nur zum Testen, ob die 13/14 Werte falsch berechnet werden (genauso wie auch hier der letzte Teil), weil beim Plotten so viel bessere Sachen mit den gleichen neuen Parametern rauskommen'''

g_AB = x0[0]
g_BA = x0[1]
alpha = 0.4  # alphaAB =alphaBA

# i = 6  #7.Wert in Tabelle, also Platz 6 in array
def rechne(Ti_mod, i):
   #hier danach jetzt einfach die gesamte Formel, aber in dem drin die Variable Ti_mod, nahc der das Nullstellenproblem gelöst werden soll --> geht nur so, wenn tau_AB, G_AB,... nicht davor
    T_new = (xB0[i] ** 2) * ((g_BA / (R * Ti_mod) * (math.exp(-alpha * g_BA / (R * Ti_mod)) / (xA0[i] + xB0[i] * math.exp(-alpha * g_BA / (R * Ti_mod)))) ** 2 + g_AB / (R * Ti_mod) * math.exp(-alpha * g_AB / (R * Ti_mod)) / (xA0[i] * math.exp(-alpha * g_AB / (R * Ti_mod)) + xB0[i]) ** 2)) - math.log(1 / xA0[i] * math.exp(-deltahLS_R / (R * Ti_mod) * (1 - Ti_mod / T0S_LS)))
    return T_new

Ti_mod_new = np.zeros(14)
T_Werte_neu = np.zeros(14)
for i in range(len(xA0)):
    Ti_mod_new[i] = T0[i]
    T_Werte_neu[i] = (opt.fsolve(rechne, Ti_mod_new[i], args=(i)))
    # FQS_new[i] = (T0[i] - T_Werte_neu[i])**2

# FQS_new_ges = sum(FQS_new)
print("xA0", xA0, "T_einWert_neu_berechnet", T_Werte_neu)
# print("Die richtige Fehlerquadratsumme (nicht normiert) beträgt:", FQS_new_ges)