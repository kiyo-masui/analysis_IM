"""Various constants
Please add a doc string if you add further constants.

::

"""

doc = ""

doc += "  pc_cm: Parsec in cm\n"
pc_cm = 3.08568025e18

doc += "  Mpc_cm: Megaparsec in cm\n"
Mpc_cm = 3.08568025e24

doc += "  Mpc_km: Megaparsec in km\n"
Mpc_km = Mpc_cm * 1.0e-5

doc += "  angstrom_cm: Angstrom in cm\n"
angstrom_cm = 1e-8

doc += "  yr_s: a year in s\n"
yr_s = 31556926.

doc += "  Myr_s: Megayear in s\n"
Myr_s = 1.e6 * yr_s

doc += "  Gyr_s: Gigayear in s\n"
Gyr_s = 1.e9 * yr_s

doc += "  amu_g: atomic mass unit in g\n"
amu_g = 1.66053886e-24

doc += "  m_p_g: mass of a proton in g\n"
m_p_g = 1.67262158e-24

doc += "  m_H_g: mass of a hydrogen atom in g\n"
m_H_g = 1.00794 * amu_g

doc += "  m_He_g: mass of a helium atom in g\n"
m_He_g = 4.002602 * amu_g

doc += "  M_sun_g: Solar mass in grams.\n"
M_sun_g = 1.98892e33

doc += "  c_light_cm_s: Speed of light in Mpc/s (from google calculator)\n"
c_light_cm_s = 29979245800.

doc += "  c_light_Mpc_s: Speed of light in Mpc/s\n"
c_light_Mpc_s = c_light_cm_s / Mpc_cm

doc += "  c_light_Mpc_Gyr: Speed of light in Mpc/Gyr \n"
c_light_Mpc_Gyr = Gyr_s * c_light_cm_s / Mpc_cm

doc += "  H100_s: 100 km s^-1 Mpc^-1 in s^-1\n"
H100_s = 100. / Mpc_km

doc += "  G_const_Mpc_Msun_s: Gravitational constant in Mpc^3 msun^-1 s^-2\n"
G_const_Mpc_Msun_s = M_sun_g * (6.673e-8) / Mpc_cm**3.

doc += "  freq_21cm_MHz: 21cm hyperfine frequency in MHz\n"
freq_21cm_MHz = 1420.40575177

doc += "  lambda_Lya_0: Central wavelength of H Lyman-alpha in Angstroms\n"
lambda_Lya_0 = 1215.67

doc += "  sigma_T_cm: Thomson cross section in cm^2\n"
sigma_T_cm = 6.6524586e-25

__doc__ += "\n".join(sorted(doc.split("\n")))
