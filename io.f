c-----------------------------------------------------------------------
c     file io.f.
c     input and output unit declarations.
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c     module declarations.
c-----------------------------------------------------------------------
      MODULE io_mod
      IMPLICIT NONE

      INTEGER :: in_unit=1
      INTEGER :: out_unit=2
      INTEGER :: bin_unit=3
      INTEGER :: term_unit=6

      INTEGER :: binr_unit=11
      INTEGER :: binc_unit=12
      INTEGER :: binj_unit=13
      INTEGER :: binx_unit=14

      INTEGER :: deltac_out_unit=21
      INTEGER :: deltac_bin_unit=22
      INTEGER :: conv_out_unit=23
      INTEGER :: conv_bin_unit=24
      INTEGER :: x01_bin_unit=25

      INTEGER :: debug_unit=99

      END MODULE io_mod
