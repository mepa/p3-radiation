m_p = 1.67262178e-24
m_e = 9.10938291e-28

m_h = m_p + m_e
m_hplu = m_p
m_hmin = m_p + 2 * m_e

m_hel = 4 * m_p + 2 * m_e
m_hep = 4 * m_p + m_e
m_hepp = 4 * m_p

m_htwo = 2 * m_p + 2 * m_e
m_htwp = 2 * m_p + m_e

m_deut = 2 * m_p + m_e
m_dplu = 2 * m_p
m_hd = 3 * m_p + 2 * m_e

m_total = m_h + m_hplu + m_hmin + m_hel + m_hep + m_hepp + m_htwo + m_htwp + m_deut + m_dplu + m_hd + m_e
m_h_nuclei = 26 * m_p

print m_total / m_h_nuclei

print '***'

print 1.0e-6 * m_total / m_h_nuclei
