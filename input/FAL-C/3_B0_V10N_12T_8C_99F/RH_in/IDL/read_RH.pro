;PRO read_RH

; ===================================================================
; === Parameters for reading Han's code output ======================
; ===================================================================

; -- Indices read from 'prd_ad.out' file. 
; This file stores the indices of the wavelengths at which the 
; radiation field is stored while performing AD PRD calculations
nwl=0
temp=0
pi=3.1415

azi=8  ; azimuthal directions

openr,1,'prd_ad.out'
repeat begin
	readf,1,temp
	nwl=nwl+1
endrep until temp eq 0
close,1

nwl=nwl-1
openr,1,'prd_ad.out'
for i=0,nwl-1 do begin
	readf,1,temp
	if(i eq 0) then wlmin=temp
	if(i eq nwl-1) then wlmax=temp
endfor
print,wlmin,wlmax
close,1

;Indices for the spectral line
line=20         ;Line to be considered
lu=8            ;Upper level
ll=0            ;Lower level

;Index for directions
imu=1           ;Index of considered direction in RH angular grid. 

; ===================================================================
; === Open output files and read output quantities ==================
; ===================================================================

f1=readinput('input.out')            ; Needed to run readspectrum
@input.common

f2=readgeometry('geometry.out')      ; Needed to run readspectrum and readatmos
@geometry.common
@files.common
nz=geometry.ndep                ;Number of heights in the atmospheric model
nmu=geometry.nrays              ;Number of ray directions
xmu=geometry.xmu                ;Ray directions (values of mu)

; -- Electrons and hydrogen atoms
metalFile  = 'metals.out'            ; Needed to run readatmos
moleculeFile = 'molecules.out'       ;   "
backgroundFile = 'background.dat'    ;   "
f3=readatmos('atmos.out')            ; Needed to run readspectrum
@opacity.common
@atmos.common
NeH=atmos.n_elec*1.d-6    ;Number density of electrons [cm^-3]
NhH=atmos.nH[*,0]*1.d-6   ;Number density of neutral hydrogen atoms in the
                          ;ground level [cm^-3]

; -- Wavelengths
f4=readspectrum('spectrum.out')
@spectrum.common
wl=spectrum.lambda(wlmin:wlmax)*10.d0   ; wavelengths [A]
wl_vac=airtovacuum(wl*0.1)*10.
nu=2.99792458d18/wl_vac                 ;Frequencies

; -- Population
atom=readatom('atom.CA.out')
pop=(*atom.n_ptr)
NlH=reform(pop[*,ll])/1.d6 ;Number density of atoms in the lower level [cm^-3]
NuH=reform(pop[*,lu])/1.d6 ;Number density of atoms in the upper level [cm^-3]

; -- Damping
f5=readDamping(atom,'damping.CA.out')
adampH=(*atom.transition[line].Adamp_ptr)  ;Damping constant

; -- Inelastic collisions
f6=readCollisions(atom,'collrate.CA.out')
CulH=(*atom.Cij_ptr)[*,lu,ll]   ;Rate of inelastic collisions

; -- Continuum
f7=openopacity('opacity.out')
print,"-- Reading background for mu=",xmu(imu)
k_cH=dblarr(nz,nwl) ;continuum total (absorption+scattering) opacity
k_sH=k_cH           ;continuum scattering opacity
eps_cH=k_cH         ;continuum thermal emission

;Read the aforementioned continuum quantities at the selected wavelengths 
;and for the selected value of mu (imu) with the routine 'readopacity'.
;Note: the dependence of the continuum on the angle is generally negligible
for i = 0, nwl-1 do begin
	j=wlmin+i
	readopacity,j,imu
	k_cH(*,i)=chi_c * 1.d-2  ;cont. total opacity [cm^-1]
	k_sH(*,i)=scatt * 1.d-2  ;cont. scatt. opacity [cm^-1]
	eps_cH(*,i)=eta_c * 10.  ;cont. therm. emiss. [erg cm^-3 s^-1 Hz^-1 sr^-1]
endfor

; -- Radiation field (J00 and J20)  
J00=readfullj(J20)
J0=dblarr(nz,nwl)
J2=J0
for i = 0, nwl-1 do begin
	  j=wlmin+i
      J0(*,i)=J00(*,j) * 1.d3  ;cgs units
	  J2(*,i)=J20(*,j) * 1.d3  ;cgs units
endfor

; -- Full radiation field I(wl,mu,ht)  
; Store radiation field in variable int(wl,mu,ht)
; Wavelengths: from shortest to longest
; Directions: from mu closest to -1 to mu closest to 1
; Heights: from top to bottom of the atmosphere
temp=dblarr(2*nmu*nwl*nz)
openr,1,'Imu.dat'
	readu,1,temp
close,1

int=dblarr(nwl,2*nmu,nz)
for i=0,nwl-1 do begin
	for mu=0,nmu-1 do begin
		offset=2*(i*nmu+mu)*nz
		int[i,nmu-mu-1,*]=temp[offset:offset+nz-1]*1.d3    ;away from obs
		int[i,nmu+mu,*]=temp[offset+nz:offset+2*nz-1]*1.d3 ;to obs
	endfor
endfor


;--- Emergent intensity (for all positive mus)
stI=spectrum.I[wlmin:wlmax,*]*1.d3  ; cgs units

; ===================================================================
; === Read elastic collisions rate from Han's code output ===========
; ===================================================================
openr,1,'qelast.out'
readf,1,levu
if(levu ne lu) then begin
	print, 'Wrong Qelast data'
	stop
endif
readf,1,levl
if(levl ne ll) then begin
	print, 'Wrong Qelast data'
	stop
endif
readf,1,grad
qelast=dblarr(nz)
readf,1,qelast
close,1


; ===================================================================
; === Read/calculate depolarizing rate D^(2) ========================
; ===================================================================
;D2=0.5*qelast
D2=0.d0*qelast

; ===================================================================
; === Open and read atmospheric model ===============================
; ===================================================================
openr,1,'../../Atmos/FALC93_70.atmos'
str=strarr(6)
readf,1,str
np=0
readf,1,np
str2=strarr(2)
readf,1,str2
temp_a=dblarr(5,np)
readf,1,temp_a
str3=strarr(3)
readf,1,str3
temp_b=dblarr(6,np)
readf,1,temp_b
close,1


height=reform(temp_a[0,*])     ;Atmospheric height [km]
tmp=reform(temp_a[1,*])        ;Temperature [K]
NeF=reform(temp_a[2,*])        ;Number density of electrons [cm^-3]
NhF=reform(temp_b[0,*])        ;Number density of HI atoms [cm^-3]
vmic=reform(temp_a[4,*])*1.d5  ;Microturbulent velocity [cm/s]
vz=reform(temp_a[3,*])         ; velocity_z [km/s]

incl_v=FLTARR(nz)               ;Inclination_v
incl_v[*]=pi
azim_v=FLTARR(nz)          ; chi_V

Bg=FLTARR(nz)              ;magnetic_field
incl_B=FLTARR(nz)          ;Inclination_B
azim_B=FLTARR(nz)          ;chi_B

; ===================================================================
; === Write input files for polarized AD PRD code ===================
; ===================================================================

openw,1,'output/dimensions.dat'
printf,1,' NZ   NTh[GL]   NCh[Tr]   NF              Note: GL=Gauss-Legendre & Tr=trapezoidal'
printf,1,nz,2*nmu,azi,nwl,format='(I3,3x,I2,7x,I2,8x,I3)'
close,1

openw,1,'output/atmosphere.dat'
;printf,1,nz,format='(i2)'
printf,1,format='("Height[km]",1x,"Temp[K]",3x,"Vmic[cm/s]",4x,"Damp",10x,"Nl[cm-3]",6x,"Nu[cm-3]",6x,"Cul[s-1]",6x,"Qel[s-1]")'
for i=0,nz-1 do begin
    printf,1,height[i],tmp[i],vmic[i],adampH[i],NlH[i],NuH[i],CulH[i],$
		Qelast[i],format='(F8.2,2x,F8.1,9(2x,E12.5))'
endfor
close,1

openw,1,'output/bulk_velocity.dat'
printf,1,format='("V[kms-1]",3x,"theta[rad] ",1x,"chi[rad]")'
for i=0,nz-1 do begin
    printf,1,vz[i],incl_v[i],azim_v[i],format='(F8.4,2x,F8.4,2x,F8.4)'
endfor
close,1


openw,1,'output/magnetic_field.dat'
printf,1,format='(2x,"B[G]",6x,"theta[rad] ",1x,"chi[rad]")'
for i=0,nz-1 do begin
    printf,1,Bg[i],incl_B[i],azim_B[i],format='(F8.4,2x,F8.4,2x,F8.4)'
endfor
close,1


openw,1,'output/continuum/continuum_tot_opac.dat'
printf,1,'cont. total opacity [cm^-1]'
for iz=0,nz-1 do begin
    printf,1,reform(k_cH[iz,*]),format='('+string(nwl)+'(2x,E12.5))'
endfor
close,1

openw,1,'output/continuum/continuum_scat_opac.dat'
printf,1,'cont. scatt. opacity [cm^-1]'
for iz=0,nz-1 do begin
	printf,1,reform(k_sH[iz,*]),format='('+string(nwl)+'(2x,E12.5))'
endfor
close,1

openw,1,'output/continuum/continuum_therm_emiss.dat'
printf,1,'cont. therm. emiss. [erg cm^-3 s^-1 Hz^-1 sr^-1]'
for iz=0,nz-1 do begin
    printf,1,reform(eps_cH[iz,*]),format='('+string(nwl)+'(2x,E12.5))'
endfor
close,1

openw,1,'output/frequency.dat'
printf,1,'Wav. (air)[A]     Freq. [s-1]'
;printf,1,'[A]           [s-1]'
for i=0,nwl-1 do begin
	printf,1,wl[i],nu[i],format='(F10.5,8x,E15.9)'
endfor
close,1


int_com=dblarr(nwl,2*nmu*azi,nz) ;defines array to include theta and chi; can certainly be improved...


for mu=0,7 do begin
   int_com[*,mu,*]=int[*,0,*]
endfor
for mu=8,15 do begin
   int_com[*,mu,*]=int[*,1,*]
endfor
for mu=16,23 do begin
   int_com[*,mu,*]=int[*,2,*]
endfor
for mu=24,31 do begin
   int_com[*,mu,*]=int[*,3,*]
endfor
for mu=32,39 do begin
   int_com[*,mu,*]=int[*,4,*]
endfor
for mu=40,47 do begin
   int_com[*,mu,*]=int[*,5,*]
endfor
for mu=48,55 do begin
   int_com[*,mu,*]=int[*,6,*]
endfor
for mu=56,63 do begin
   int_com[*,mu,*]=int[*,7,*]
endfor
for mu=64,71 do begin
   int_com[*,mu,*]=int[*,8,*]
endfor
for mu=72,79 do begin
   int_com[*,mu,*]=int[*,9,*]
endfor
for mu=80,87 do begin
   int_com[*,mu,*]=int[*,10,*]
endfor
for mu=88,95 do begin
   int_com[*,mu,*]=int[*,11,*]
endfor



for iz=0,nz-1 do begin
filename=string(iz,format='("output/radiation/StokesI_0",i2,".res")')
if(iz lt 10) then filename=string(iz,format='("output/radiation/StokesI_00",i1,".res")')
openw,1,filename
	for mu=0,(2*nmu)*azi-1 do begin
		printf,1,reform(int_com[*,mu,iz-1]),format='('+string(nwl)+'(2x,E12.5))'
	endfor
close,1
endfor

; Needed to initialize AA PRD code
openw,1,'output/Jrad.dat'
printf,1,nz,nwl,format='(i2,2x,i3)'
for j=0,nz-1 do begin
	printf,1,reform(J0[j,*]),format='('+string(nwl)+'(2x,E12.5))'
endfor
for j=0,nz-1 do begin
	printf,1,reform(J2[j,*]),format='('+string(nwl)+'(2x,E12.5))'
endfor
close,1

; Emergent radiation calculated by RH (for reference)
openw,1,'output/emergent.out'
printf,1,nwl,nmu,format='(i3,2x,i2)'
	for mu=0,nmu-1 do begin
		printf,1,reform(stI[*,mu]),format='('+string(nwl)+'(2x,E12.5))'
	endfor
close,1

stop

END
