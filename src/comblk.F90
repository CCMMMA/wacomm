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


      integer nsources, i_source(600), j_source(600), k_source(600),  &
                     id_source(600), npartsperhour(600), mode(600),   &
                     source_start(600), source_end(600)


      end module
