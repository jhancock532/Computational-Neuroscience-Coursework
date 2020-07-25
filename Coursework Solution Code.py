import random as rnd
import matplotlib.pyplot as plt
import numpy as np

# PART A

# Question One
#Simulate an integrate and fire model with the following parameters for 1 s:

#timeConstant = 10 * 10 ** -3 # ms #membrane time constant        Tm
#leakVoltage = -70 * 10 ** -3 #mV #leak potential                 El
#restVoltage = -70 * 10 ** -3 #mV #rest voltage / reset voltage   Vr
#thresholdVoltage = -40 * 10 ** -3 # mV #voltage threshold        Vth
#membraneResistance = 10 * 10 ** 6 #M #membrane resistance.       Rm
#externalCurrent = 3.1 * 10 ** -9 #nA. #external current signal   Ie
#
#timestep = 0.25 * 10 ** -3 #ms
#
#### Verbose Method ###
#
##derivativeAtZero = (externalCurrent * membraneResistance) / timeConstant
##currentDerivative = derivativeAtZero
##
##for i in range(1, 4000):
##    #calculate V given derivative and x
##    voltage = currentDerivative * timestep + VValues[i-1]
##    
##    #if the voltage is above the threshold, reduce it.
##    if voltage > thresholdVoltage:
##        voltage = restVoltage
##        
##    VValues.append(voltage)
##
##    #find the next derivative
##    currentDerivative = ((restVoltage - voltage) + externalCurrent * membraneResistance) / timeConstant
#
#### Simplified Method ###
#
#dt = timestep
#El = leakVoltage
#IeRm = externalCurrent * membraneResistance
#Tm = timeConstant
#Vrest = restVoltage
#Vth = thresholdVoltage
#
#tValues = np.arange(0, 1, timestep)
#VValues = [restVoltage]
#
#dV = ((El - VValues[0]) + IeRm) * (dt / Tm)
#
#for i in range(1, 4000):
#    vPrev = VValues[i-1]
#    
#    v = vPrev + dV
#    
#    if v > Vth:
#        v = Vrest
#    
#    dV = ((El - v) + IeRm) * (dt / Tm)
#    
#    VValues.append(v)
#
### Plotting ###
#plt.rcParams.update({'font.size': 16})
#fig, ax = plt.subplots(figsize=(15,9))
#ax.plot(tValues, VValues, 'k', linewidth=1)
#ax.set(xlabel='Time (s)', ylabel='Voltage (V)',
#       title='Integrate-and-Fire Simulation')
#ax.grid()
#fig.savefig("question one A figure.png")
#plt.show()

## Question Two
##Neuron Parameters
#Tm = 20 * 10 ** -3
#El = -70 * 10 ** -3
#Vrest = -80 * 10 ** -3
#Vth = -54 * 10 ** -3
#RmIe = 18 * 10 ** -3
##Synapse Parameters
#Rmgs = 0.15
#deltaS = 0.5 #a spike causes s to increase instantly by this amount
#Ts = 10 * 10 ** -3
#
##Simulate two cases:
##Case One
#Es = 0 #synapses are excitatory
##Case Two
##Es = -80 * 10 ** -3 #synapses are inhibitory
#
##Initial neuron voltages selected randomly
#v1 = np.random.uniform(low=Vrest, high=Vth, size=(1,))[0]
#v2 = np.random.uniform(low=Vrest, high=Vth, size=(1,))[0]
#
#timestep = 0.25 * 10 ** -3 #ms
#tValues = np.arange(0, 1, timestep)
#dt = timestep
#
## Set initial voltages
#v1Values = [v1]
#v2Values = [v2]
#s1Values = [0]
#s2Values = [0]
#
#n1Spiked = False
#n2Spiked = False
#
#for i in range(1, 4000):
#    v1Prev = v1Values[i-1]
#    v2Prev = v2Values[i-1]
#    sv1Prev = s1Values[i-1]
#    sv2Prev = s2Values[i-1]
#    
#    if n2Spiked:
#        sv1 = sv1Prev - sv1Prev * (dt / Ts) + deltaS
#        n2Spiked = False
#    else:
#        sv1 = sv1Prev - sv1Prev * (dt / Ts)
#        
#    if n1Spiked:
#        sv2 = sv2Prev - sv2Prev * (dt / Ts) + deltaS
#        n1Spiked = False
#    else:
#        sv2 = sv2Prev - sv2Prev * (dt / Ts)
#    
#    RmIs1 = Rmgs * (Es - v1Prev) * sv1 #The voltage of synapse one.
#    RmIs2 = Rmgs * (Es - v2Prev) * sv2 #The voltage of synapse two.
#    
#    dV1 = ((El - v1Prev) + RmIe + RmIs1) * (dt / Tm)
#    dV2 = ((El - v2Prev) + RmIe + RmIs2) * (dt / Tm)
#    
#    v1 = v1Prev + dV1
#    v2 = v2Prev + dV2
#    
#    if v1 > Vth:
#        v1 = Vrest
#        n1Spiked = True
#    
#    if v2 > Vth:
#        v2 = Vrest
#        n2Spiked = True
#        
#    s1Values.append(sv1)
#    s2Values.append(sv2)
#    v1Values.append(v1)
#    v2Values.append(v2)
#
#plt.rcParams.update({'font.size': 16})
#fig, ax = plt.subplots(figsize=(15,9))
#ax.plot(tValues, v1Values, 'k', linewidth=1)
#ax.plot(tValues, v2Values, 'b', linewidth=1)
##ax.plot(tValues, s1Values, 'k', linewidth=1)
##ax.plot(tValues, s2Values, 'b', linewidth=1)
#ax.set(xlabel='Time (s)', ylabel='Voltage (V)',
#       title='Neurons with Excitatory Synapses')
#ax.grid()
#fig.savefig("Task A Question Two.png")
#plt.show()

### Part B
# Question 1

#dt = 0.25 * 10 ** -3
#El = -65 * 10 ** -3
#Tm = 10 * 10 ** -3
#Rm = 100 * 10 ** 6
#Vrest = -65 * 10 ** -3
#Vth = -50 * 10 ** -3
#
#numSynapses = 40
#Ts = 2 * 10 ** -3
#Es = 0
#deltaS = 0.5
#firingRate = 15
#
##The initial peak conductance
#initialgBar = 4 * 10 ** -9 #nanoSiemens (A value of conductance)
#
#gBarArray = [initialgBar] * 40
#sArray = []
#for x in range(0, 40):
#    sArray.append([0])
#RmIsArray = [0] * 40
#
#tValues = np.arange(0, 1, dt)
#nVValues = [Vrest]
#
#numSpikes = 0
#for i in range(1, 4000):
#    synapseVoltageTally = 0
#    
#    nVPrev = nVValues[i-1]
#    
#    #Test the state of all the synapses.
#    #Check if they have fired then update their s term
#    for synapseNum in range(numSynapses):
#        
#        sPrev = sArray[synapseNum][i-1]
#        
#        if np.random.uniform() < firingRate * dt:
#            sV = sPrev - sPrev * (dt / Ts) + deltaS
#        else:
#            sV = sPrev - sPrev * (dt / Ts)
#            
#        sArray[synapseNum].append(sV)
#
#        synapseVoltageTally += (gBarArray[synapseNum] * sV)
#    
#    RmIs = Rm * synapseVoltageTally * (Es - nVPrev)
#
#    dnV = ((El - nVPrev) + RmIs) * (dt / Tm)
#    
#    nV = nVPrev + dnV
#    
#    if nV > Vth:
#        nV = Vrest
#        numSpikes += 1
#        
#    nVValues.append(nV)
#
#titleText = "Neuron with 40 Synapses (" + str(numSpikes) + " Spikes)"
#
#plt.rcParams.update({'font.size': 16})
#fig, ax = plt.subplots(figsize=(15,9))
#ax.plot(tValues, nVValues, 'k', linewidth=1)
#
#ax.set(xlabel='Time (s)', ylabel='Voltage (V)',
#       title=titleText)
#ax.grid()
#plt.show()


### Part B 
### Question Two & Three

#dt = 0.25 * 10 ** -3
#El = -65 * 10 ** -3
#Tm = 10 * 10 ** -3
#Rm = 100 * 10 ** 6
#Vrest = -65 * 10 ** -3
#Vth = -50 * 10 ** -3
#
#numSynapses = 40
#Ts = 2 * 10 ** -3
#Es = 0
#deltaS = 0.5
#firingRate = 10
#
##The initial peak conductance
#initialgBar = 4 * 10 ** -9 #nanoSiemens (A value of conductance)
#tpreArray = [0] * 40 #Double check on Reddit
#sampleGBar = []
#tpost = -1000 #0
#
#Aplus = 0.2 * 10 ** -9
#Aminus = 0.25 * 10 ** -9
#Tplus = 20 * 10 ** -3
#Tminus = 20 * 10 ** -3
#
#simulationLength = 300 #in seconds.
#tValues = np.arange(0, simulationLength, dt)
#
##for w in range(10, 22, 2):
#w = 20
#    
#print("\n\nNEW RUN RESULTS FOR " + str(w) + "HZ\n")
#firingRate = w
#
#simulationFrames = int(simulationLength * (1 / dt))
#
#gBarArray = [initialgBar] * 40
#sVArray = []
#for x in range(0, 40):
#    sVArray.append([0])
#
#nVValues = [Vrest]
#
#STDP_ON_FLAG = True
#POST_SPIKE_OCCURED = False
#
#numSpikes = 0
#spikeBinResults = []
#
#for i in range(1, simulationFrames+1):
#    
#    if (dt * i) % 10 == 0:
#        spikeBinResults.append(numSpikes / 10)
#        numSpikes = 0
#    
#    synapseVoltageTally = 0
#    
#    nVPrev = nVValues[i-1]
#    
#    #Test the state of all the synapses.
#    #Check if they have fired then update their s term
#    for synapseNum in range(numSynapses):
#        PRE_SPIKE_OCCURED = False
#        sVPrev = sVArray[synapseNum][i-1]
#        
#        if np.random.uniform() < firingRate * dt:
#            PRE_SPIKE_OCCURED = True
#            sV = sVPrev - sVPrev * (dt / Ts) + deltaS
#        else:
#            sV = sVPrev - sVPrev * (dt / Ts)
#            
#        sVArray[synapseNum].append(sV)
#        
#        if STDP_ON_FLAG:
#            # Update the gBar weights
#            gBarPrev = gBarArray[synapseNum]
#            
#            tDelta = tpost - tpreArray[synapseNum]
#            
#            gBar = gBarPrev
#            
#            if PRE_SPIKE_OCCURED:
#                tpreArray[synapseNum] = dt * i
#                gBar = gBarPrev - Aminus * np.exp(-abs(tDelta)/Tminus)
#                PRE_SPIKE_OCCURED = False
#                
#            if POST_SPIKE_OCCURED:
#                gBar = gBarPrev + Aplus * np.exp(-abs(tDelta)/Tplus)
#            
#            if gBar < 0:
#                gBar = 0
#            if gBar > initialgBar:
#                gBar = initialgBar
#                
#            gBarArray[synapseNum] = gBar
#            
#            if synapseNum == 5:
#                sampleGBar.append(gBar)
#
#        synapseVoltageTally += (gBarArray[synapseNum] * sV)
#        
#    POST_SPIKE_OCCURED = False
#    
#    RmIs = Rm * synapseVoltageTally * (Es - nVPrev)
#
#    dnV = ((El - nVPrev) + RmIs) * (dt / Tm)
#    
#    nV = nVPrev + dnV
#    
#    if nV > Vth:
#        nV = Vrest
#        tpost = dt * i
#        POST_SPIKE_OCCURED = True
#        numSpikes += 1
#        
#    nVValues.append(nV)
#
#titleText = "Neuron Simulation: "
#if STDP_ON_FLAG:
#    titleText += "STDP Flag On"
#else:
#    titleText += "STDP Flag Off"
#
#print("Logging of Synapse Weights")
#for i in range(0, 40):
#    print(str(gBarArray[i]))
#    
#print("Spike Bin Results")
#for i in range(0, len(spikeBinResults)):
#    print(spikeBinResults[i])
#
#plt.rcParams.update({'font.size': 16})
#fig, ax = plt.subplots(figsize=(15,9))
#ax.plot(tValues, sampleGBar, 'k', linewidth=1) #nVValues
#
#ax.set(xlabel='Time (s)', ylabel='Voltage (V)',
#       title=titleText)
#ax.grid()
#plt.show()

##### PART B ######
## Question Four ##

dt = 0.25 * 10 ** -3
El = -65 * 10 ** -3
Tm = 10 * 10 ** -3
Rm = 100 * 10 ** 6
Vrest = -65 * 10 ** -3
Vth = -50 * 10 ** -3

numSynapses = 40
Ts = 2 * 10 ** -3
Es = 0
deltaS = 0.5
firingRate = 10

#The initial peak conductance
initialgBar = 4 * 10 ** -9 #nanoSiemens (A value of conductance)
tpreArray = [0] * 40
sampleGBar = []
tpost = -1000 #0

Aplus = 0.2 * 10 ** -9
Aminus = 0.25 * 10 ** -9
Tplus = 20 * 10 ** -3
Tminus = 20 * 10 ** -3

simulationLength = 300 #in seconds.
tValues = np.arange(0, simulationLength, dt)

B = 20
averageFiringRate = 20

simulationFrames = int(simulationLength * (1 / dt))

gBarArray = [initialgBar] * 40
sVArray = []
for x in range(0, 40):
    sVArray.append([0])

nVValues = [Vrest]

STDP_ON_FLAG = True
POST_SPIKE_OCCURED = False

numSpikes = 0
spikeBinResults = []

for i in range(1, simulationFrames+1):
    
    if (dt * i) % 10 == 0:
        spikeBinResults.append(numSpikes / 10)
        numSpikes = 0
    
    synapseVoltageTally = 0
    
    nVPrev = nVValues[i-1]
    
    #Test the state of all the synapses.
    #Check if they have fired then update their s term
    for synapseNum in range(numSynapses):
        PRE_SPIKE_OCCURED = False
        sVPrev = sVArray[synapseNum][i-1]
        
        if np.random.uniform() < (averageFiringRate + B * np.sin(2 * np.pi * 10 * i * dt)) * dt:
            PRE_SPIKE_OCCURED = True
            sV = sVPrev - sVPrev * (dt / Ts) + deltaS
        else:
            sV = sVPrev - sVPrev * (dt / Ts)
            
        sVArray[synapseNum].append(sV)
        
        if STDP_ON_FLAG:
            # Update the gBar weights
            gBarPrev = gBarArray[synapseNum]
            
            tDelta = tpost - tpreArray[synapseNum]
            
            gBar = gBarPrev
            
            if PRE_SPIKE_OCCURED:
                tpreArray[synapseNum] = dt * i
                gBar = gBarPrev - Aminus * np.exp(-abs(tDelta)/Tminus)
                PRE_SPIKE_OCCURED = False
                
            if POST_SPIKE_OCCURED:
                gBar = gBarPrev + Aplus * np.exp(-abs(tDelta)/Tplus)
            
            if gBar < 0:
                gBar = 0
            if gBar > initialgBar:
                gBar = initialgBar
                
            gBarArray[synapseNum] = gBar
            
            if synapseNum == 5:
                sampleGBar.append(gBar)

        synapseVoltageTally += (gBarArray[synapseNum] * sV)
        
    POST_SPIKE_OCCURED = False
    
    RmIs = Rm * synapseVoltageTally * (Es - nVPrev)

    dnV = ((El - nVPrev) + RmIs) * (dt / Tm)
    
    nV = nVPrev + dnV
    
    if nV > Vth:
        nV = Vrest
        tpost = dt * i
        POST_SPIKE_OCCURED = True
        numSpikes += 1
        
    nVValues.append(nV)

titleText = "Neuron Simulation: "
if STDP_ON_FLAG:
    titleText += "STDP Flag On"
else:
    titleText += "STDP Flag Off"

print("Logging of Synapse Weights")
for i in range(0, 40):
    print(str(gBarArray[i]))
    
print("Spike Bin Results")
for i in range(0, len(spikeBinResults)):
    print(spikeBinResults[i])

plt.rcParams.update({'font.size': 16})
fig, ax = plt.subplots(figsize=(15,9))
ax.plot(tValues, nVValues, 'k', linewidth=1)

ax.set(xlabel='Time (s)', ylabel='Voltage (V)',
       title=titleText)
ax.grid()
plt.show()

