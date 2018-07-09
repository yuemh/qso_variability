# The Continuum Variability of SDSS-RM Quasars

## Intro (Previous Studies)
There has been models, codes, observations on quasar continuum variability. The most important observation is MacLeod et al. (2010), who use SDSS Stripe 82 data, with longest baseline ~9 year.

As far as I know, MacLeod et al. (2010) work is the basic of most of the knowledge of this topic. But I believe that we can dig further.

## Aim

Investigate the following topics:
- What are the factors that influences the DRW parameters? Or, put it in another way, what can the measured DRW parameters tell us about quasars?
- Can DRW model describe the short timescale variation of quasars? What can the short timescale variations tell us about the accretion disk
- Summarize the distribution of quasar variabilities

## Sample

Quasars: SDSS_RM quasars (849  quasars)

Photometric data should include:

- SDSS photometry
- SDSS DR7 spectrum
- CFHT legacy survey (g band)
- Pan-Starrs (shall we wait for DR2?)
- SDSS DR12 spectrum
- SDSS-RM
- BASS (??)

ALL quasars should have ~13 year time baseline, with number of photometric points >80.

Need to check with Ian and Xiaohui about the origin of data.

### Combining different surveys
- SDSS photometry: PSF magnitude
- SDSS DR7 spectrum: integrate the spectrum and do a typical aperture correction, assuming a typical seeing.
- CFHT legacy survey (g band): use SDSS-RM to correct Ian's data
- Pan-Starrs (shall we wait for DR2?): Aperture mag. (Pan-Starrs has performed aperture correction)
- SDSS DR12 spectrum: integrate the spectrum and do a typical aperture correction, assuming a typical seeing.
- SDSS-RM: used as the standard
- BASS: use SDSS-RM to correct Ian's data

## Analyze

Use JAVELIN to measure the parameters.

First analyze the influence of rest-frame wavelength. Compare with rpevious studies. (Consider what does this mean).

Note the influence of bias. especially on tau.
