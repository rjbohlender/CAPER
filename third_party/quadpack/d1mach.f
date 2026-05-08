      double precision function d1mach(i)
      integer i
      if (i .eq. 1) then
        d1mach = tiny(1.0d0)
      else if (i .eq. 2) then
        d1mach = huge(1.0d0)
      else if (i .eq. 3) then
        d1mach = epsilon(1.0d0) / 2.0d0
      else if (i .eq. 4) then
        d1mach = epsilon(1.0d0)
      else if (i .eq. 5) then
        d1mach = log10(2.0d0)
      else
        d1mach = 0.0d0
      endif
      return
      end
