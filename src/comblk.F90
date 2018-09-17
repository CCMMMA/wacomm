!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      module comblk

      integer totpart

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      integer idead,iouter,ip,                              &
              isorgente,jsorgente,ksorgente,ifluxsorgente


      real deltax,deltay
      integer npart
!--------------------------------------------------------------------------


      real xpart, ypart, zpart, health
      real health0, tau0, survprob
      integer pstatus

      allocatable xpart(:), ypart(:), zpart(:), health(:), tpart(:)
      allocatable pstatus(:)


      integer nsources, i_source(100), j_source(100), k_source(100),  &
                     id_source(100), npartsperhour(100), mode(100),   &
                     source_start(100), source_end(100)


      end module
