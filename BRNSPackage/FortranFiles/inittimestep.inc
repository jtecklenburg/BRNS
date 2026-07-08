c    $Id: inittimestep.inc 16 2026-07-08 10:56:47Z tecklenburg $
      BLOCK DATA InitTimeStep
        include 'timestep.inc'
        DATA dtold/1.d-6/, ttol/5.d-2/, tstep/1.0d0/, maxconc/0.d0/
      END