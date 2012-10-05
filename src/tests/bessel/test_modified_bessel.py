# This script provides a function rat(x) that gives rational approximation to
# the function besseli(4+0.5, x) / exp(x) for any x > 0. It plots an error plot
# showing, that both absolute and relative errors are less than 1e-15 against
# the exact answer calculated using mpmath.  The rational approximation was
# calculated using Mathematica, a sample script is given in the comment below.

from numpy import linspace
from math import sinh, cosh, exp, sqrt, pi
from sympy.mpmath import besseli

from common import make_plots

def Ik3(x):
    if x < 0.4:
        r = x**4/105 + x**6/1890 + x**8/83160 + x**10/6486480 + \
            x**12/778377600 + x**14/132324192e3
        r = r / exp(x)
    else:
        r = -(15/x**3 + 6/x)*(sinh(x)/exp(x)) + (15/x**2 + 1)*(cosh(x)/exp(x))
    r = r * sqrt(2/(pi*x))
    return r

def Ik4(x):
    if x < 0.2:
        r = x**5/945 + x**7/20790 + x**9/1081080 + x**11/97297200 + \
            x**13/132324192e2
        r = r/exp(x)
    elif x < 1.7:
        r =  (5.1962956007054264229527112912e-17 -  \
                1.7091983797257865101834560981e-15*x +  \
                2.54701673734823215056998758576e-14*x**2 -  \
                2.28250221664547749050764003136e-13*x**3 +  \
                1.37741730213761458713328353583e-12*x**4 +  \
                0.00105820105225973559701867830214*x**5 +  \
                0.00061204753548652964602354522028*x**6 +  \
                0.0000122798136500667925710228809483*x**7 +  \
                0.0000103551941240484592046206146607*x**8)/ \
            (1 + 0.578384903100157953572414461702*x -  \
                0.0338500781802291050601812074314*x**2 -  \
                0.0165046449338415288740437462715*x**3 +  \
                0.000664629591462297572628738277521*x**4 +  \
                0.000244501438485835213668147591168*x**5 -  \
                0.0000102233345389484523128091507614*x**6 -  \
                2.37680351471245713258958333627e-6*x**7 +  \
                1.7722924118248819944246867464e-7*x**8)
        r = r/exp(x)
    elif x < 4:
        r = (-1.53382275599031353712671171402e-7 + \
                1.02445790975414106564424700453e-6*x - \
                3.20332766010782394861350740562e-6*x**2 + \
                6.22984887800335387305228161855e-6*x**3 - \
                8.44422432663036167543192411535e-6*x**4 + \
                0.00106667316609330486086491285357*x**5 + \
                0.000296879521931548405158141021032*x**6 + \
                9.16131185419642061414089694738e-6*x**7 + \
                4.84398673121860851791923295866e-6*x**8)/\
            (1 + 0.286715277505499007952029608584*x - \
                0.0405239647887260390437917379491*x**2 - \
                0.00665840458518766453585847552188*x**3 + \
                0.000270261982447901375808140464959*x**4 + \
                0.000271851192449937947349585685247*x**5 - \
                0.0000431403377516821653299354827079*x**6 + \
                2.65420161107604023280035849057e-6*x**7 - \
                6.20361275935376828687251028769e-8*x**8)
        r = r/exp(x)
    elif x < 10:
        # Produced by:
        # FortranForm[MiniMaxApproximation[((105/x^4 + 45/x^2 +
        #   1)*Sinh[x]-(105/x^3 + 10/x)*Cosh[x])/Exp[x], {x, {4, 10}, 8, 8},
        #   WorkingPrecision->30]]
        r = (0.000395502959013236968661582656143 - \
                0.001434648369704841686633794071*x + \
                0.00248783474583503473135143644434*x**2 - \
                0.00274477921388295929464613063609*x**3 + \
                0.00216275018107657273725589740499*x**4 - \
                0.000236779926184242197820134964535*x**5 + \
                0.0000882030507076791807159699814428*x**6 - \
                4.62078105288798755556136693122e-6*x**7 + \
                8.23671374777791529292655504214e-7*x**8)/\
            (1 + 0.504839286873735708062045336271*x + \
                0.176683950009401712892997268723*x**2 + \
                0.0438594911840609324095487447279*x**3 + \
                0.00829753062428409331123592322788*x**4 + \
                0.00111693697900468156881720995034*x**5 + \
                0.000174719963536517752971223459247*x**6 + \
                7.22885338737473776714257581233e-6*x**7 + \
                1.64737453771748367647332279826e-6*x**8)
    elif x < 20:
        r = (105/x**4 + 45/x**2 + 1)*(sinh(x)/exp(x)) - \
                (105/x**3 + 10/x)*(cosh(x)/exp(x))
    else:
        r = (105/x**4 - 105/x**3 + 45/x**2 - 10/x + 1)/2

    r = r * sqrt(2/(pi*x))
    return r

xx = linspace(1e-10, 40, 10000)
yf = [(besseli(4+0.5, x) / exp(x)) for x in xx]
yrat = [Ik4(x) for x in xx]

make_plots(xx, yf, yrat)
