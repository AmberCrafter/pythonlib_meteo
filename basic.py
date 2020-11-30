# =================================================================== #
# Author: AmberCraft
# License: MIT
# Version: 1.0.0
# Type: toolbox
# Publish date: 2020-09-06
# =================================================================== #
# <ChangeLog>
# <Version: 1.0.0>
# This is the first public version, maybe has some bug. If yuo find any
# bug, please let me known to fix it.
# =================================================================== #
# Declear Parameter/Variable
# variableName: type <unit> [option] -- description and note
# ------------------------------------------------------------------- #
# @basic variable
# longitude: float <deg>
# latitude: float <deg>
# altitude: float <deg>
# temperature: float <degC>
# dTemperature: float <degC> -- dry temperature
# wTemperature: float <degC> -- wet temperature
# rHumidity: float <%> -- relative humidity
# pressure: float <hPa>
# windSpeed: float <m/s>
# windDir: float <deg>
# rainfall: float <mm>
# ------------------------------------------------------------------- #
# @derive variable
# eTemperature: float <K> -- equivalent temperature
# epTemperature: float <K> -- equivalent potential temperature
# dewTemperature: float <degC> -- dewpoint temperature
# mixRatio: float <1> -- mixing ratio
# pTemperature: float <K> -- potential temperature
# spHumidity: float <g/kg> -- specific humidity
# vapor: float <hPa> -- water vapor pressure
# vTemperature: float <K> -- virtural temperature
# vpTemperature: float <K> -- virtural potential temperature
# ------------------------------------------------------------------- #
# @other parameter
# AGL: float <m> -- above ground level height
# density_air <kg/m3>
# spLatent: float <kJ/kg> -- specific latent heat
# triTemperature_water: float <K> -- water triple point temperature
# =================================================================== #
# Declear function
# functionName -- discription and not
# ------------------------------------------------------------------- #
# dewpoint_temperature
# mixing_ratio
# equivalent_potential_temperature
# exp_air_density
# exp_gravity
# equivalent_temperature
# exp_specific_latent_evap
# moist_static_energy
# potential_temperature
# saturation_vapor_pressure
# specific_humidity
# vapor_pressure
# virtural_potential_temperature
# virtural_temperature
# =================================================================== #
# Function description formation
'''
# Input Parameter
@type variable: type <unit> [option] -- discription and note
# Option Input Parameter
@type variable -> {
    @type variable: type <unit> [option] -- discription and not
}
# Output Parameter
@type variable: type <unit> [option] -- discription and not
@result

# Function
# Discription
# Ref and Src
'''
# @type = [
#     basi = basic variable, 
#     deri = derive variable,
#     disc = discription parameter,
#     para = other parameter
# ]
# =================================================================== #

def _ndarray_float(x):
    import numpy as np
    return np.array(x,dtype=float)

def dewpoint_temperature(temperautre,rHumidity):
    '''
    # Input Parameter:
    @basi temperature: float <degC>
    @basi rHumidity: float <%>

    # Output Parameter:
    @deri dewTemperature: float <degC>

    # Discription
    1. effective range
        0℃ < T < 60℃
        1% < RH < 100%
        0℃ < T_{d} < 50℃
    2. depend on numpy to calculate array

    # Src: wiki
    '''
    import numpy as np
    temperautre=_ndarray_float(temperautre)
    rHumidity=_ndarray_float(rHumidity)

    a=17.27
    b=237.7
    gamma=a*temperautre/(b+temperautre)+np.log(rHumidity/100)
    return b*gamma/(a-gamma)

def saturation_vapor_pressure(temperature,phase='water'):
    '''
    # Input Parameter:
    @basi temperature: float <degC>
    @disc phase: str <1> ['water','ice'] -- defined what the water phase at 0℃

    # Output Parameter:
    @deri satVapor: float <hPa> -- saturation vapor pressure

    # Discription
    1. effective range
        -80℃ < T < 50℃
    2. depend on numpy to calculate array

    # Formula
        Arden Buck equations -> es = 6.1121*exp((18.678-(temperature/234.5))*(temperature/(257.14+temperature)))
        Goff Gratch equations ->     # Not implement
            phase='water' -> es = 10**(
                C1*(1-triTemperature_water/(temperature+273.15)) +
                C2*log((temperature+273.15)/triTemperature_water) +
                C3*10**(-4)*(1-10**(C4*(temperature/triTemperature_water)**(-1))) +
                C5*10**(-3)*(10**(C6*(1-triTemperature_water/temperature))-1) +
                C7)
            C1 = 10.79586
            C2 = -5.02808
            C3 = 1.50474
            C4 = -8.20602
            C5 = 0.42873
            C6 = 4.76955
            C7 = 0.7861183
            phase='ice' -> es = 10**(
                C1*(triTemperature_water/temperature-1) +
                C2*log(triTemperature_water/temperature) +
                C3(1-temperature/triTemperature_water) +
                C4
            )
            C1 = -9.906936
            C2 = -3.56654
            C3 = 0.876817
            C4 = 0.7861183
        triTemperature_water = 273.16 <K>
    # Ref: 
        Arden Buck equations
        Goff Gratch equation(1946)
    '''
    import numpy as np
    def liquid(temperature):
        return 6.1121*np.exp((18.678-(temperature/234.5))*(temperature/(257.14+temperature)))
    def solid(temperature):
        # incorrect
        return 6.1121*np.exp((18.678-(temperature/234.5))*(temperature/(257.14+temperature)))
    
    phase_list={
        'solid':0,
        's':0,
        'ice':0,

        'liquid':1,
        'l':1,
        'water':1
    }
    
    temperature=_ndarray_float(temperature)

    if phase_list[phase]==1:
        Tl=np.array(temperature,dtype=float)
        Tl[Tl<0]=np.nan
        Ts=np.array(temperature,dtype=float)
        Ts[Ts>=0]=np.nan
    else:
        Tl=np.array(temperature,dtype=float)
        Tl[Tl<=0]=np.nan
        Ts=np.array(temperature,dtype=float)
        Ts[Ts>0]=np.nan
    Tl=liquid(Tl); Ts=solid(Ts)
    return np.nansum([Tl,Ts],axis=0)

def vapor_pressure(temperature=None,rHumidity=None,method='rh',dTemperature=None,wTemperature=None,Pressure=None,*args,**keywords):
    '''
    # Input Parameter
    @basi temperature: float <degC>
    @basi rHumidity: float <%> -- relative humidity
    @disc method: str <1> [default='rh', 'dw'] -- defined the calculation method, dw mean dry and wet bulb temperature measurement method
        {
            if method=='dw' -> keywords:{
                @basi dTemperature: float <degC>
                @basi wTemperature: float <degC>
                @basi pressure: float <hPa>
            }
        }
    @disc phase: str <1> [default='water','ice'] -- defined what the water phase at 0℃

    # Output Parameter
    @deri vapor pressure: float <hPa>

    # Formula
    method='rh' ->
        e=es*rh
    method='wd' ->
        phase='water' -> e=es-0.5*(dTemperature-wTemperature)*Pressure/1013.25  (es is water saturation pressure)
        phase='ice'   -> e=es-0.44*(dTemperature-wTemperature)*Pressure/1013.25 (es is ice saturation pressure)
    method='Ferrel'  # Not implement
        phase='water' -> e=es+C1*Pressure*(dTemperature-wTemperature)*(1+C2*(wTemperature+273.15))
            C1 = -6.6*10**(-4)
            C2 = 0.00115
        phase='ice'   -> e=es+C1*Pressure*(dTemperature-wTemperature)*(1+C2*(wTemperature+273.15))
            C1 = -5.973*10**(-4)
            C2 = 0.00115

    # Src: https://zhidao.baidu.com/question/373004649658108244.html
    '''
    if method in ['rh']:
        return saturation_vapor_pressure(temperature,*args,**keywords)*rHumidity/100
    if method in ['dw','wd']:
        if 'phase' in keywords.keys():
            if keywords['phase']=='water': return saturation_vapor_pressure(dtemperature,*args,**keywords)-0.5*(dTemperature-wTemperature)*Pressure/1013.25
            if keywords['phase']=='ice': return saturation_vapor_pressure(dtemperature,*args,**keywords)-0.44*(dTemperature-wTemperature)*Pressure/1013.25
        else:
            return saturation_vapor_pressure(dTemperature,*args,**keywords)-0.5*(dTemperature-wTemperature)*Pressure/1013.25
            

def mixing_ratio(pressure,vapor=None,*args,**keywords):
    '''
    # Input Parameter:
    @basi pressure: float <hPa>
    @deri vapor: float <hPa>

    #Option Input Parameter
    @deri vapor->keywords:{
        @basi temperature: float <degC>
        @basi rHumidity: float <%>
        @basi phase: str [default='water']
    }
    
    # Output Parameter:
    @deri mixRatio: float <1> -- mixing ratio

    # Function
    mixRatio = 0.622 * (vapor/(pressure-vapor))

    # Src: AMS Glossary
    '''
    if vapor==None: vapor=vapor_pressure(**keywords)
    return 0.622*(vapor/(pressure-vapor))

def specific_humidity(pressure,vapor=None,*args,**keywords):
    '''
    # Input Parameter:
    @basi pressure: float <hPa>
    @deri vapor: float <hPa>
    
    # Option Input Parameter
    @deri vapor -> {
        @basi temperature: float
        @basi rHumidity: float <%>
        @disc phase: str <1> [defaul='water','ice'] -- defined what the water phase at 0℃
    }
    
    # output:
    @deri specific humidity: float <g/kg>
    
    depend on numpy to calculate array
    src: wiki
    '''
    # return 0.622*(vapor/pressure)
    if vapor==None: vapor=vapor_pressure(*args,**keywords)
    q = mixing_ratio(pressure,vapor)
    return q/(1+q)*1000

def potential_temperature(temperature,pressure):
    '''
    # Input Parameter
    @basi temperature: float <K>
    @basi pressure: float <hPa>

    # Output Parameter
    @deri pTemperature: float <K> -- potential temperature

    #Formula
    theta = T*(P0/P)**(R/Cp)
    P0=1000 <hPa>
    '''
    R =0.287    #unit: kJ/(kg*K)
    Cp=1.005    #unit: kJ/(kg*K)
    P0=1000     #unit: hPa
    return temperature*(P0/pressure)**(R/Cp)

def virtural_temperature(temperature,pressure,vapor=None,*args,**keywords):
    '''
    # Input Parameter
    @basi temperature: float <degC>
    @basi pressure: float  <hPa>
    @deri vapor: float <hPa> -- vapor pressure

    # Option Input Parameter
    @deri vapor -> {
        *@basi temperature: float <degC>
        @basi rHumidity: float <%>
        @disc phase: str <1> [defaul='water','ice'] -- defined what the water phase at 0℃
    }

    # Ouput Parameter
    @deri vTemperature: float <K> -- virture temperature

    #Formula
    Tv = T/(1-(e/p)*(1-epsilon))
    e: vapor pressure
    epsilon=Rd/Rv=Mv/Md~=0.622

    # Ref: https://en.wikipedia.org/wiki/Virtual_temperature
    '''
    epsilon=0.622
    if vapor==None: vapor=vapor_pressure(temperature,**keywords)
    temperature = temperature+273.15
    return temperature/(1-(vapor/pressure)*(1-epsilon))

def virtural_potential_temperature(*args,**keywords):
    '''
    # Input Parameter
    @basi temperature: float <degC>
    @basi pressure: float <hPa>
    @deri vapor: float <hPa> -- vapor pressure

    # Option Input Parameter
    @deri vapor -> {
        @basi temperature: float <degC>
        @basi rHumidity: float <%>
    }

    # Ouput Parameter
    @deri vpTemperature: float <K> -- virture potential temperature

    #Formula
    theta = Tv*(P0/P)**(R/Cp)
    Tv = T/(1-(e/p)*(1-epsilon))
    e: vapor pressure
    epsilon=Rd/Rv=Mv/Md~=0.622

    # Ref: https://en.wikipedia.org/wiki/Virtual_temperature
    '''
    return potential_temperature(temperature=virtural_temperature(**keywords),pressure=keywords['pressure'])

def exp_specific_latent_evap(temperature,phase='water'):
    '''
    # Input Parameter
    @basi temperature: float <degC>
    @desc phase: str <1> [default='water','ice']

    #Output Parameter
    @para spLatent: float <kJ/kg> -- specific latent heat

    # Ref: https://en.wikipedia.org/wiki/Latent_heat
    '''
    if phase=='water':
        return (2500.8-2.36*temperature+0.0016*temperature**2-0.00006*temperature**3)
    if phase=='ice':
        return (2831.1-0.29*temperature-0.004*temperature**2)

def exp_gravity(latitude,altitude):
    '''
    # Input Parameter
    @basi latitude: float <deg>
    @basi altitude: float <m>

    # Ouput Parameter
    @para gravity: float <m/s2>

    # Formula
    g ~= g0 * (1 + 0.0052884*(sin(latitude)**2) - 0.0000059*(sin(2*latitude)**2)) - 0.000003086*altitude
    g0 ~= 9.78046 <m/s2>

    # Discription
    1. Depend on numpy

    # Ref
    https://zh.wikipedia.org/wiki/%E9%87%8D%E5%8A%9B%E5%8A%A0%E9%80%9F%E5%BA%A6
    '''
    import numpy as np
    latitude=np.deg2rad(latitude)
    g0 = 9.78046
    return g0 * (1 + 0.0052884*(np.sin(latitude)**2) - 0.0000059*(np.sin(2*latitude)**2)) - 0.000003086*altitude

def exp_air_density(temperautre,pressure):
    '''
    # Input Parameter
    @basi temperature: float <degC>
    @basi pressure: float <hPa>

    # Output Parameter
    @para density_air: float <kg/m3>

    # Formula
    density_air = 1.293 * (pressure/1033.6) * (273.15/(273.15+temperature))
    '''
    return 1.293*(pressure/1033.6)*(273.15/(273.15+temperautre))

def moist_static_energy(temperature,pressure,AGL,mixRatio=None,*args,**keywords):
    '''
    # Input Parameter
    @basi temperature: float <degC>
    @basi pressure: float <hPa>
    @para AGL: float <m> -- above ground level height
    @deri mixRatio: float <1>

    # Option Input Parameter
    @deri mixRatio -> {
        @basi temperature: float <degC>
        @basi rHumidity: float <%>
        @basi pressure: float <hPa>
    }

    # Ouput Parameter
    @result moist static energy: float <kJ/kg>

    # Formula
    S = Cp*T + g*z + Lv*q
    z: above ground level heigh (AGL) <m>
    q: spHumidity <g/kg>
    Lv: latent heat of vaporation <kJ/g>

    # Ref
    https://en.wikipedia.org/wiki/Moist_static_energy
    '''
    rho = lambda Tair, Press: 1.293*(Press/1033.6)*(273.15/(273.15+Tair))
    Cp=1.005 #unit: kJ/(kg*K)
    Lv=exp_specific_latent_evap(temperature)/1000
    g=exp_gravity(24.968498,132)
    if mixRatio==None: mixRatio=mixing_ratio(temperature=temperature,pressure=pressure,*args,**keywords)
    return Cp*temperature + g*AGL/rho(temperature,pressure) + Lv*mixRatio

def equivalent_temperature(temperature,rHumidity,pressure):
    '''
    # Input Parameter
    @basi temperature: float <degC>
    @basi rHumidity: float <%>
    @basi pressure: float <hPa>

    # Ouput Parameter
    @deri relative temperature: float <K>

    # Formula
    Te = T*exp((Lv*mixingRatio)/(Cp*Td))
    Lv: latent heat of vaporation <kJ/g>
    Td: dewpoint temperature

    # Discription
    1. Depend on numpy
    
    # Ref
    https://zh.wikipedia.org/wiki/%E7%9B%B8%E7%95%B6%E4%BD%8D%E6%BA%AB
    '''
    import numpy as np
    Cp=1.005 #unit: kJ/(kg*K)
    Lv=exp_specific_latent_evap(temperature)/1000
    dewTemperature = dewpoint_temperature(temperature,rHumidity)
    mixRatio=mixing_ratio(temperature=temperature,rHumidity=rHumidity,pressure=pressure)
    temperature+=273.15
    return temperature*np.exp((Lv*mixRatio)/(Cp*dewTemperature))

def equivalent_potential_temperature(*args,**keywords):
    '''
    # Input Parameter
    @basi temperature: float <degC>
    @basi rHumidity: float <%>
    @basi pressure: float <hPa>

    # Ouput Parameter
    @deri epTemperature: float <K> -- equivalent potential temperature

    # Formula
    theta_e = Te*(P0/P)**(R/Cp)

    # Ref
    https://zh.wikipedia.org/wiki/%E7%9B%B8%E7%95%B6%E4%BD%8D%E6%BA%AB
    '''
    return potential_temperature(temperature=equivalent_temperature(**keywords),pressure=keywords['pressure'])

if __name__ == "__main__":
    Temperature = 25
    

    satPressure = saturation_vapor_pressure(25)
    vapPressure = vapor_pressure(dTemperature=25, wTemperature=20, method='dw', Pressure=1000)

    print('stop')
