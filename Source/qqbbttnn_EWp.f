c ============================= qqb_bbbtatann_EWp.f ============================
c file originally generated by madgraph
c edited by: stefano moretti and declan millar
c d.millar@soton.ac.uk
c ============================= sqqb_bbbtatann_EWp =============================
c returns amplitude squared summed/avg over colors and helicities
c for the point in phase space p1,p2,p3,p4,p5,p6,p7,p8  
c for process : q q -> b b l+ l- v_l v_l via (via A ,Z ,{Z'})
c =============================== start  =======================================
      real*8 function sqqb_bbbtatann_EWp(iq,gZpq,gZpt,Zpmass,Zpwidth
     &                                         ,p1,p2,p3,p4,p5,p6,p7,p8)
c note iq (initial quark) has now been added to the arguments of this program.
c 1 = down, 2 = up, 3 = strange, 4 = charm, 5 = bottom, 6 = top
c ----------------------------- declarations -----------------------------------
      implicit none

c constants 
      integer nexternal ! number of external legs   
      integer ncomb ! number of helicity combinations
      parameter ( nexternal=8, ncomb=256 )

c arguments
      integer iq ! incoming quark type (up/down)
      real*8 gZpq(2,5),gZpt(2,5) ! Zp couplings to quark (up/down) and top 
      real*8 Zpmass(5), Zpwidth(5) ! Zp mass and width
      real*8 p1(0:3),p2(0:3)
     & ,p3(0:3),p4(0:3),p5(0:3),p6(0:3),p7(0:3),p8(0:3)
  
c local variables 
      integer nhel(nexternal,ncomb),ntry
      real*8 t
      real*8 qqb_bbbtatann_EWp
      integer ihel ! helicity corresponding to p_i
      logical goodhel(ncomb)
      data goodhel/ncomb*.false./
      data ntry/0/
c all possible helicity combinations
      data (nhel(ihel,  1),ihel=1,8) / -1, -1, -1, -1, -1, -1, -1, -1/
      data (nhel(ihel,  2),ihel=1,8) / -1, -1, -1, -1, -1, -1, -1,  1/
      data (nhel(ihel,  3),ihel=1,8) / -1, -1, -1, -1, -1, -1,  1, -1/
      data (nhel(ihel,  4),ihel=1,8) / -1, -1, -1, -1, -1, -1,  1,  1/
      data (nhel(ihel,  5),ihel=1,8) / -1, -1, -1, -1, -1,  1, -1, -1/
      data (nhel(ihel,  6),ihel=1,8) / -1, -1, -1, -1, -1,  1, -1,  1/
      data (nhel(ihel,  7),ihel=1,8) / -1, -1, -1, -1, -1,  1,  1, -1/
      data (nhel(ihel,  8),ihel=1,8) / -1, -1, -1, -1, -1,  1,  1,  1/
      data (nhel(ihel,  9),ihel=1,8) / -1, -1, -1, -1,  1, -1, -1, -1/
      data (nhel(ihel, 10),ihel=1,8) / -1, -1, -1, -1,  1, -1, -1,  1/
      data (nhel(ihel, 11),ihel=1,8) / -1, -1, -1, -1,  1, -1,  1, -1/
      data (nhel(ihel, 12),ihel=1,8) / -1, -1, -1, -1,  1, -1,  1,  1/
      data (nhel(ihel, 13),ihel=1,8) / -1, -1, -1, -1,  1,  1, -1, -1/
      data (nhel(ihel, 14),ihel=1,8) / -1, -1, -1, -1,  1,  1, -1,  1/
      data (nhel(ihel, 15),ihel=1,8) / -1, -1, -1, -1,  1,  1,  1, -1/
      data (nhel(ihel, 16),ihel=1,8) / -1, -1, -1, -1,  1,  1,  1,  1/
      data (nhel(ihel, 17),ihel=1,8) / -1, -1, -1,  1, -1, -1, -1, -1/
      data (nhel(ihel, 18),ihel=1,8) / -1, -1, -1,  1, -1, -1, -1,  1/
      data (nhel(ihel, 19),ihel=1,8) / -1, -1, -1,  1, -1, -1,  1, -1/
      data (nhel(ihel, 20),ihel=1,8) / -1, -1, -1,  1, -1, -1,  1,  1/
      data (nhel(ihel, 21),ihel=1,8) / -1, -1, -1,  1, -1,  1, -1, -1/
      data (nhel(ihel, 22),ihel=1,8) / -1, -1, -1,  1, -1,  1, -1,  1/
      data (nhel(ihel, 23),ihel=1,8) / -1, -1, -1,  1, -1,  1,  1, -1/
      data (nhel(ihel, 24),ihel=1,8) / -1, -1, -1,  1, -1,  1,  1,  1/
      data (nhel(ihel, 25),ihel=1,8) / -1, -1, -1,  1,  1, -1, -1, -1/
      data (nhel(ihel, 26),ihel=1,8) / -1, -1, -1,  1,  1, -1, -1,  1/
      data (nhel(ihel, 27),ihel=1,8) / -1, -1, -1,  1,  1, -1,  1, -1/
      data (nhel(ihel, 28),ihel=1,8) / -1, -1, -1,  1,  1, -1,  1,  1/
      data (nhel(ihel, 29),ihel=1,8) / -1, -1, -1,  1,  1,  1, -1, -1/
      data (nhel(ihel, 30),ihel=1,8) / -1, -1, -1,  1,  1,  1, -1,  1/
      data (nhel(ihel, 31),ihel=1,8) / -1, -1, -1,  1,  1,  1,  1, -1/
      data (nhel(ihel, 32),ihel=1,8) / -1, -1, -1,  1,  1,  1,  1,  1/
      data (nhel(ihel, 33),ihel=1,8) / -1, -1,  1, -1, -1, -1, -1, -1/
      data (nhel(ihel, 34),ihel=1,8) / -1, -1,  1, -1, -1, -1, -1,  1/
      data (nhel(ihel, 35),ihel=1,8) / -1, -1,  1, -1, -1, -1,  1, -1/
      data (nhel(ihel, 36),ihel=1,8) / -1, -1,  1, -1, -1, -1,  1,  1/
      data (nhel(ihel, 37),ihel=1,8) / -1, -1,  1, -1, -1,  1, -1, -1/
      data (nhel(ihel, 38),ihel=1,8) / -1, -1,  1, -1, -1,  1, -1,  1/
      data (nhel(ihel, 39),ihel=1,8) / -1, -1,  1, -1, -1,  1,  1, -1/
      data (nhel(ihel, 40),ihel=1,8) / -1, -1,  1, -1, -1,  1,  1,  1/
      data (nhel(ihel, 41),ihel=1,8) / -1, -1,  1, -1,  1, -1, -1, -1/
      data (nhel(ihel, 42),ihel=1,8) / -1, -1,  1, -1,  1, -1, -1,  1/
      data (nhel(ihel, 43),ihel=1,8) / -1, -1,  1, -1,  1, -1,  1, -1/
      data (nhel(ihel, 44),ihel=1,8) / -1, -1,  1, -1,  1, -1,  1,  1/
      data (nhel(ihel, 45),ihel=1,8) / -1, -1,  1, -1,  1,  1, -1, -1/
      data (nhel(ihel, 46),ihel=1,8) / -1, -1,  1, -1,  1,  1, -1,  1/
      data (nhel(ihel, 47),ihel=1,8) / -1, -1,  1, -1,  1,  1,  1, -1/
      data (nhel(ihel, 48),ihel=1,8) / -1, -1,  1, -1,  1,  1,  1,  1/
      data (nhel(ihel, 49),ihel=1,8) / -1, -1,  1,  1, -1, -1, -1, -1/
      data (nhel(ihel, 50),ihel=1,8) / -1, -1,  1,  1, -1, -1, -1,  1/
      data (nhel(ihel, 51),ihel=1,8) / -1, -1,  1,  1, -1, -1,  1, -1/
      data (nhel(ihel, 52),ihel=1,8) / -1, -1,  1,  1, -1, -1,  1,  1/
      data (nhel(ihel, 53),ihel=1,8) / -1, -1,  1,  1, -1,  1, -1, -1/
      data (nhel(ihel, 54),ihel=1,8) / -1, -1,  1,  1, -1,  1, -1,  1/
      data (nhel(ihel, 55),ihel=1,8) / -1, -1,  1,  1, -1,  1,  1, -1/
      data (nhel(ihel, 56),ihel=1,8) / -1, -1,  1,  1, -1,  1,  1,  1/
      data (nhel(ihel, 57),ihel=1,8) / -1, -1,  1,  1,  1, -1, -1, -1/
      data (nhel(ihel, 58),ihel=1,8) / -1, -1,  1,  1,  1, -1, -1,  1/
      data (nhel(ihel, 59),ihel=1,8) / -1, -1,  1,  1,  1, -1,  1, -1/
      data (nhel(ihel, 60),ihel=1,8) / -1, -1,  1,  1,  1, -1,  1,  1/
      data (nhel(ihel, 61),ihel=1,8) / -1, -1,  1,  1,  1,  1, -1, -1/
      data (nhel(ihel, 62),ihel=1,8) / -1, -1,  1,  1,  1,  1, -1,  1/
      data (nhel(ihel, 63),ihel=1,8) / -1, -1,  1,  1,  1,  1,  1, -1/
      data (nhel(ihel, 64),ihel=1,8) / -1, -1,  1,  1,  1,  1,  1,  1/
      data (nhel(ihel, 65),ihel=1,8) / -1,  1, -1, -1, -1, -1, -1, -1/
      data (nhel(ihel, 66),ihel=1,8) / -1,  1, -1, -1, -1, -1, -1,  1/
      data (nhel(ihel, 67),ihel=1,8) / -1,  1, -1, -1, -1, -1,  1, -1/
      data (nhel(ihel, 68),ihel=1,8) / -1,  1, -1, -1, -1, -1,  1,  1/
      data (nhel(ihel, 69),ihel=1,8) / -1,  1, -1, -1, -1,  1, -1, -1/
      data (nhel(ihel, 70),ihel=1,8) / -1,  1, -1, -1, -1,  1, -1,  1/
      data (nhel(ihel, 71),ihel=1,8) / -1,  1, -1, -1, -1,  1,  1, -1/
      data (nhel(ihel, 72),ihel=1,8) / -1,  1, -1, -1, -1,  1,  1,  1/
      data (nhel(ihel, 73),ihel=1,8) / -1,  1, -1, -1,  1, -1, -1, -1/
      data (nhel(ihel, 74),ihel=1,8) / -1,  1, -1, -1,  1, -1, -1,  1/
      data (nhel(ihel, 75),ihel=1,8) / -1,  1, -1, -1,  1, -1,  1, -1/
      data (nhel(ihel, 76),ihel=1,8) / -1,  1, -1, -1,  1, -1,  1,  1/
      data (nhel(ihel, 77),ihel=1,8) / -1,  1, -1, -1,  1,  1, -1, -1/
      data (nhel(ihel, 78),ihel=1,8) / -1,  1, -1, -1,  1,  1, -1,  1/
      data (nhel(ihel, 79),ihel=1,8) / -1,  1, -1, -1,  1,  1,  1, -1/
      data (nhel(ihel, 80),ihel=1,8) / -1,  1, -1, -1,  1,  1,  1,  1/
      data (nhel(ihel, 81),ihel=1,8) / -1,  1, -1,  1, -1, -1, -1, -1/
      data (nhel(ihel, 82),ihel=1,8) / -1,  1, -1,  1, -1, -1, -1,  1/
      data (nhel(ihel, 83),ihel=1,8) / -1,  1, -1,  1, -1, -1,  1, -1/
      data (nhel(ihel, 84),ihel=1,8) / -1,  1, -1,  1, -1, -1,  1,  1/
      data (nhel(ihel, 85),ihel=1,8) / -1,  1, -1,  1, -1,  1, -1, -1/
      data (nhel(ihel, 86),ihel=1,8) / -1,  1, -1,  1, -1,  1, -1,  1/
      data (nhel(ihel, 87),ihel=1,8) / -1,  1, -1,  1, -1,  1,  1, -1/
      data (nhel(ihel, 88),ihel=1,8) / -1,  1, -1,  1, -1,  1,  1,  1/
      data (nhel(ihel, 89),ihel=1,8) / -1,  1, -1,  1,  1, -1, -1, -1/
      data (nhel(ihel, 90),ihel=1,8) / -1,  1, -1,  1,  1, -1, -1,  1/
      data (nhel(ihel, 91),ihel=1,8) / -1,  1, -1,  1,  1, -1,  1, -1/
      data (nhel(ihel, 92),ihel=1,8) / -1,  1, -1,  1,  1, -1,  1,  1/
      data (nhel(ihel, 93),ihel=1,8) / -1,  1, -1,  1,  1,  1, -1, -1/
      data (nhel(ihel, 94),ihel=1,8) / -1,  1, -1,  1,  1,  1, -1,  1/
      data (nhel(ihel, 95),ihel=1,8) / -1,  1, -1,  1,  1,  1,  1, -1/
      data (nhel(ihel, 96),ihel=1,8) / -1,  1, -1,  1,  1,  1,  1,  1/
      data (nhel(ihel, 97),ihel=1,8) / -1,  1,  1, -1, -1, -1, -1, -1/
      data (nhel(ihel, 98),ihel=1,8) / -1,  1,  1, -1, -1, -1, -1,  1/
      data (nhel(ihel, 99),ihel=1,8) / -1,  1,  1, -1, -1, -1,  1, -1/
      data (nhel(ihel,100),ihel=1,8) / -1,  1,  1, -1, -1, -1,  1,  1/
      data (nhel(ihel,101),ihel=1,8) / -1,  1,  1, -1, -1,  1, -1, -1/
      data (nhel(ihel,102),ihel=1,8) / -1,  1,  1, -1, -1,  1, -1,  1/
      data (nhel(ihel,103),ihel=1,8) / -1,  1,  1, -1, -1,  1,  1, -1/
      data (nhel(ihel,104),ihel=1,8) / -1,  1,  1, -1, -1,  1,  1,  1/
      data (nhel(ihel,105),ihel=1,8) / -1,  1,  1, -1,  1, -1, -1, -1/
      data (nhel(ihel,106),ihel=1,8) / -1,  1,  1, -1,  1, -1, -1,  1/
      data (nhel(ihel,107),ihel=1,8) / -1,  1,  1, -1,  1, -1,  1, -1/
      data (nhel(ihel,108),ihel=1,8) / -1,  1,  1, -1,  1, -1,  1,  1/
      data (nhel(ihel,109),ihel=1,8) / -1,  1,  1, -1,  1,  1, -1, -1/
      data (nhel(ihel,110),ihel=1,8) / -1,  1,  1, -1,  1,  1, -1,  1/
      data (nhel(ihel,111),ihel=1,8) / -1,  1,  1, -1,  1,  1,  1, -1/
      data (nhel(ihel,112),ihel=1,8) / -1,  1,  1, -1,  1,  1,  1,  1/
      data (nhel(ihel,113),ihel=1,8) / -1,  1,  1,  1, -1, -1, -1, -1/
      data (nhel(ihel,114),ihel=1,8) / -1,  1,  1,  1, -1, -1, -1,  1/
      data (nhel(ihel,115),ihel=1,8) / -1,  1,  1,  1, -1, -1,  1, -1/
      data (nhel(ihel,116),ihel=1,8) / -1,  1,  1,  1, -1, -1,  1,  1/
      data (nhel(ihel,117),ihel=1,8) / -1,  1,  1,  1, -1,  1, -1, -1/
      data (nhel(ihel,118),ihel=1,8) / -1,  1,  1,  1, -1,  1, -1,  1/
      data (nhel(ihel,119),ihel=1,8) / -1,  1,  1,  1, -1,  1,  1, -1/
      data (nhel(ihel,120),ihel=1,8) / -1,  1,  1,  1, -1,  1,  1,  1/
      data (nhel(ihel,121),ihel=1,8) / -1,  1,  1,  1,  1, -1, -1, -1/
      data (nhel(ihel,122),ihel=1,8) / -1,  1,  1,  1,  1, -1, -1,  1/
      data (nhel(ihel,123),ihel=1,8) / -1,  1,  1,  1,  1, -1,  1, -1/
      data (nhel(ihel,124),ihel=1,8) / -1,  1,  1,  1,  1, -1,  1,  1/
      data (nhel(ihel,125),ihel=1,8) / -1,  1,  1,  1,  1,  1, -1, -1/
      data (nhel(ihel,126),ihel=1,8) / -1,  1,  1,  1,  1,  1, -1,  1/
      data (nhel(ihel,127),ihel=1,8) / -1,  1,  1,  1,  1,  1,  1, -1/
      data (nhel(ihel,128),ihel=1,8) / -1,  1,  1,  1,  1,  1,  1,  1/
      data (nhel(ihel,129),ihel=1,8) /  1, -1, -1, -1, -1, -1, -1, -1/
      data (nhel(ihel,130),ihel=1,8) /  1, -1, -1, -1, -1, -1, -1,  1/
      data (nhel(ihel,131),ihel=1,8) /  1, -1, -1, -1, -1, -1,  1, -1/
      data (nhel(ihel,132),ihel=1,8) /  1, -1, -1, -1, -1, -1,  1,  1/
      data (nhel(ihel,133),ihel=1,8) /  1, -1, -1, -1, -1,  1, -1, -1/
      data (nhel(ihel,134),ihel=1,8) /  1, -1, -1, -1, -1,  1, -1,  1/
      data (nhel(ihel,135),ihel=1,8) /  1, -1, -1, -1, -1,  1,  1, -1/
      data (nhel(ihel,136),ihel=1,8) /  1, -1, -1, -1, -1,  1,  1,  1/
      data (nhel(ihel,137),ihel=1,8) /  1, -1, -1, -1,  1, -1, -1, -1/
      data (nhel(ihel,138),ihel=1,8) /  1, -1, -1, -1,  1, -1, -1,  1/
      data (nhel(ihel,139),ihel=1,8) /  1, -1, -1, -1,  1, -1,  1, -1/
      data (nhel(ihel,140),ihel=1,8) /  1, -1, -1, -1,  1, -1,  1,  1/
      data (nhel(ihel,141),ihel=1,8) /  1, -1, -1, -1,  1,  1, -1, -1/
      data (nhel(ihel,142),ihel=1,8) /  1, -1, -1, -1,  1,  1, -1,  1/
      data (nhel(ihel,143),ihel=1,8) /  1, -1, -1, -1,  1,  1,  1, -1/
      data (nhel(ihel,144),ihel=1,8) /  1, -1, -1, -1,  1,  1,  1,  1/
      data (nhel(ihel,145),ihel=1,8) /  1, -1, -1,  1, -1, -1, -1, -1/
      data (nhel(ihel,146),ihel=1,8) /  1, -1, -1,  1, -1, -1, -1,  1/
      data (nhel(ihel,147),ihel=1,8) /  1, -1, -1,  1, -1, -1,  1, -1/
      data (nhel(ihel,148),ihel=1,8) /  1, -1, -1,  1, -1, -1,  1,  1/
      data (nhel(ihel,149),ihel=1,8) /  1, -1, -1,  1, -1,  1, -1, -1/
      data (nhel(ihel,150),ihel=1,8) /  1, -1, -1,  1, -1,  1, -1,  1/
      data (nhel(ihel,151),ihel=1,8) /  1, -1, -1,  1, -1,  1,  1, -1/
      data (nhel(ihel,152),ihel=1,8) /  1, -1, -1,  1, -1,  1,  1,  1/
      data (nhel(ihel,153),ihel=1,8) /  1, -1, -1,  1,  1, -1, -1, -1/
      data (nhel(ihel,154),ihel=1,8) /  1, -1, -1,  1,  1, -1, -1,  1/
      data (nhel(ihel,155),ihel=1,8) /  1, -1, -1,  1,  1, -1,  1, -1/
      data (nhel(ihel,156),ihel=1,8) /  1, -1, -1,  1,  1, -1,  1,  1/
      data (nhel(ihel,157),ihel=1,8) /  1, -1, -1,  1,  1,  1, -1, -1/
      data (nhel(ihel,158),ihel=1,8) /  1, -1, -1,  1,  1,  1, -1,  1/
      data (nhel(ihel,159),ihel=1,8) /  1, -1, -1,  1,  1,  1,  1, -1/
      data (nhel(ihel,160),ihel=1,8) /  1, -1, -1,  1,  1,  1,  1,  1/
      data (nhel(ihel,161),ihel=1,8) /  1, -1,  1, -1, -1, -1, -1, -1/
      data (nhel(ihel,162),ihel=1,8) /  1, -1,  1, -1, -1, -1, -1,  1/
      data (nhel(ihel,163),ihel=1,8) /  1, -1,  1, -1, -1, -1,  1, -1/
      data (nhel(ihel,164),ihel=1,8) /  1, -1,  1, -1, -1, -1,  1,  1/
      data (nhel(ihel,165),ihel=1,8) /  1, -1,  1, -1, -1,  1, -1, -1/
      data (nhel(ihel,166),ihel=1,8) /  1, -1,  1, -1, -1,  1, -1,  1/
      data (nhel(ihel,167),ihel=1,8) /  1, -1,  1, -1, -1,  1,  1, -1/
      data (nhel(ihel,168),ihel=1,8) /  1, -1,  1, -1, -1,  1,  1,  1/
      data (nhel(ihel,169),ihel=1,8) /  1, -1,  1, -1,  1, -1, -1, -1/
      data (nhel(ihel,170),ihel=1,8) /  1, -1,  1, -1,  1, -1, -1,  1/
      data (nhel(ihel,171),ihel=1,8) /  1, -1,  1, -1,  1, -1,  1, -1/
      data (nhel(ihel,172),ihel=1,8) /  1, -1,  1, -1,  1, -1,  1,  1/
      data (nhel(ihel,173),ihel=1,8) /  1, -1,  1, -1,  1,  1, -1, -1/
      data (nhel(ihel,174),ihel=1,8) /  1, -1,  1, -1,  1,  1, -1,  1/
      data (nhel(ihel,175),ihel=1,8) /  1, -1,  1, -1,  1,  1,  1, -1/
      data (nhel(ihel,176),ihel=1,8) /  1, -1,  1, -1,  1,  1,  1,  1/
      data (nhel(ihel,177),ihel=1,8) /  1, -1,  1,  1, -1, -1, -1, -1/
      data (nhel(ihel,178),ihel=1,8) /  1, -1,  1,  1, -1, -1, -1,  1/
      data (nhel(ihel,179),ihel=1,8) /  1, -1,  1,  1, -1, -1,  1, -1/
      data (nhel(ihel,180),ihel=1,8) /  1, -1,  1,  1, -1, -1,  1,  1/
      data (nhel(ihel,181),ihel=1,8) /  1, -1,  1,  1, -1,  1, -1, -1/
      data (nhel(ihel,182),ihel=1,8) /  1, -1,  1,  1, -1,  1, -1,  1/
      data (nhel(ihel,183),ihel=1,8) /  1, -1,  1,  1, -1,  1,  1, -1/
      data (nhel(ihel,184),ihel=1,8) /  1, -1,  1,  1, -1,  1,  1,  1/
      data (nhel(ihel,185),ihel=1,8) /  1, -1,  1,  1,  1, -1, -1, -1/
      data (nhel(ihel,186),ihel=1,8) /  1, -1,  1,  1,  1, -1, -1,  1/
      data (nhel(ihel,187),ihel=1,8) /  1, -1,  1,  1,  1, -1,  1, -1/
      data (nhel(ihel,188),ihel=1,8) /  1, -1,  1,  1,  1, -1,  1,  1/
      data (nhel(ihel,189),ihel=1,8) /  1, -1,  1,  1,  1,  1, -1, -1/
      data (nhel(ihel,190),ihel=1,8) /  1, -1,  1,  1,  1,  1, -1,  1/
      data (nhel(ihel,191),ihel=1,8) /  1, -1,  1,  1,  1,  1,  1, -1/
      data (nhel(ihel,192),ihel=1,8) /  1, -1,  1,  1,  1,  1,  1,  1/
      data (nhel(ihel,193),ihel=1,8) /  1,  1, -1, -1, -1, -1, -1, -1/
      data (nhel(ihel,194),ihel=1,8) /  1,  1, -1, -1, -1, -1, -1,  1/
      data (nhel(ihel,195),ihel=1,8) /  1,  1, -1, -1, -1, -1,  1, -1/
      data (nhel(ihel,196),ihel=1,8) /  1,  1, -1, -1, -1, -1,  1,  1/
      data (nhel(ihel,197),ihel=1,8) /  1,  1, -1, -1, -1,  1, -1, -1/
      data (nhel(ihel,198),ihel=1,8) /  1,  1, -1, -1, -1,  1, -1,  1/
      data (nhel(ihel,199),ihel=1,8) /  1,  1, -1, -1, -1,  1,  1, -1/
      data (nhel(ihel,200),ihel=1,8) /  1,  1, -1, -1, -1,  1,  1,  1/
      data (nhel(ihel,201),ihel=1,8) /  1,  1, -1, -1,  1, -1, -1, -1/
      data (nhel(ihel,202),ihel=1,8) /  1,  1, -1, -1,  1, -1, -1,  1/
      data (nhel(ihel,203),ihel=1,8) /  1,  1, -1, -1,  1, -1,  1, -1/
      data (nhel(ihel,204),ihel=1,8) /  1,  1, -1, -1,  1, -1,  1,  1/
      data (nhel(ihel,205),ihel=1,8) /  1,  1, -1, -1,  1,  1, -1, -1/
      data (nhel(ihel,206),ihel=1,8) /  1,  1, -1, -1,  1,  1, -1,  1/
      data (nhel(ihel,207),ihel=1,8) /  1,  1, -1, -1,  1,  1,  1, -1/
      data (nhel(ihel,208),ihel=1,8) /  1,  1, -1, -1,  1,  1,  1,  1/
      data (nhel(ihel,209),ihel=1,8) /  1,  1, -1,  1, -1, -1, -1, -1/
      data (nhel(ihel,210),ihel=1,8) /  1,  1, -1,  1, -1, -1, -1,  1/
      data (nhel(ihel,211),ihel=1,8) /  1,  1, -1,  1, -1, -1,  1, -1/
      data (nhel(ihel,212),ihel=1,8) /  1,  1, -1,  1, -1, -1,  1,  1/
      data (nhel(ihel,213),ihel=1,8) /  1,  1, -1,  1, -1,  1, -1, -1/
      data (nhel(ihel,214),ihel=1,8) /  1,  1, -1,  1, -1,  1, -1,  1/
      data (nhel(ihel,215),ihel=1,8) /  1,  1, -1,  1, -1,  1,  1, -1/
      data (nhel(ihel,216),ihel=1,8) /  1,  1, -1,  1, -1,  1,  1,  1/
      data (nhel(ihel,217),ihel=1,8) /  1,  1, -1,  1,  1, -1, -1, -1/
      data (nhel(ihel,218),ihel=1,8) /  1,  1, -1,  1,  1, -1, -1,  1/
      data (nhel(ihel,219),ihel=1,8) /  1,  1, -1,  1,  1, -1,  1, -1/
      data (nhel(ihel,220),ihel=1,8) /  1,  1, -1,  1,  1, -1,  1,  1/
      data (nhel(ihel,221),ihel=1,8) /  1,  1, -1,  1,  1,  1, -1, -1/
      data (nhel(ihel,222),ihel=1,8) /  1,  1, -1,  1,  1,  1, -1,  1/
      data (nhel(ihel,223),ihel=1,8) /  1,  1, -1,  1,  1,  1,  1, -1/
      data (nhel(ihel,224),ihel=1,8) /  1,  1, -1,  1,  1,  1,  1,  1/
      data (nhel(ihel,225),ihel=1,8) /  1,  1,  1, -1, -1, -1, -1, -1/
      data (nhel(ihel,226),ihel=1,8) /  1,  1,  1, -1, -1, -1, -1,  1/
      data (nhel(ihel,227),ihel=1,8) /  1,  1,  1, -1, -1, -1,  1, -1/
      data (nhel(ihel,228),ihel=1,8) /  1,  1,  1, -1, -1, -1,  1,  1/
      data (nhel(ihel,229),ihel=1,8) /  1,  1,  1, -1, -1,  1, -1, -1/
      data (nhel(ihel,230),ihel=1,8) /  1,  1,  1, -1, -1,  1, -1,  1/
      data (nhel(ihel,231),ihel=1,8) /  1,  1,  1, -1, -1,  1,  1, -1/
      data (nhel(ihel,232),ihel=1,8) /  1,  1,  1, -1, -1,  1,  1,  1/
      data (nhel(ihel,233),ihel=1,8) /  1,  1,  1, -1,  1, -1, -1, -1/
      data (nhel(ihel,234),ihel=1,8) /  1,  1,  1, -1,  1, -1, -1,  1/
      data (nhel(ihel,235),ihel=1,8) /  1,  1,  1, -1,  1, -1,  1, -1/
      data (nhel(ihel,236),ihel=1,8) /  1,  1,  1, -1,  1, -1,  1,  1/
      data (nhel(ihel,237),ihel=1,8) /  1,  1,  1, -1,  1,  1, -1, -1/
      data (nhel(ihel,238),ihel=1,8) /  1,  1,  1, -1,  1,  1, -1,  1/
      data (nhel(ihel,239),ihel=1,8) /  1,  1,  1, -1,  1,  1,  1, -1/
      data (nhel(ihel,240),ihel=1,8) /  1,  1,  1, -1,  1,  1,  1,  1/
      data (nhel(ihel,241),ihel=1,8) /  1,  1,  1,  1, -1, -1, -1, -1/
      data (nhel(ihel,242),ihel=1,8) /  1,  1,  1,  1, -1, -1, -1,  1/
      data (nhel(ihel,243),ihel=1,8) /  1,  1,  1,  1, -1, -1,  1, -1/
      data (nhel(ihel,244),ihel=1,8) /  1,  1,  1,  1, -1, -1,  1,  1/
      data (nhel(ihel,245),ihel=1,8) /  1,  1,  1,  1, -1,  1, -1, -1/
      data (nhel(ihel,246),ihel=1,8) /  1,  1,  1,  1, -1,  1, -1,  1/
      data (nhel(ihel,247),ihel=1,8) /  1,  1,  1,  1, -1,  1,  1, -1/
      data (nhel(ihel,248),ihel=1,8) /  1,  1,  1,  1, -1,  1,  1,  1/
      data (nhel(ihel,249),ihel=1,8) /  1,  1,  1,  1,  1, -1, -1, -1/
      data (nhel(ihel,250),ihel=1,8) /  1,  1,  1,  1,  1, -1, -1,  1/
      data (nhel(ihel,251),ihel=1,8) /  1,  1,  1,  1,  1, -1,  1, -1/
      data (nhel(ihel,252),ihel=1,8) /  1,  1,  1,  1,  1, -1,  1,  1/
      data (nhel(ihel,253),ihel=1,8) /  1,  1,  1,  1,  1,  1, -1, -1/
      data (nhel(ihel,254),ihel=1,8) /  1,  1,  1,  1,  1,  1, -1,  1/
      data (nhel(ihel,255),ihel=1,8) /  1,  1,  1,  1,  1,  1,  1, -1/
      data (nhel(ihel,256),ihel=1,8) /  1,  1,  1,  1,  1,  1,  1,  1/
c ---------------------------------- method ------------------------------------
      sqqb_bbbtatann_EWp = 0d0
      ntry=ntry+1
      do ihel=1,ncomb
          if (goodhel(ihel) .or. ntry .lt. 10) then
             t=qqb_bbbtatann_EWp(iq,gZpq,gZpt
     &        ,Zpmass,Zpwidth,p1,p2,p3,p4,p5,p6,p7,p8,nhel(1,ihel)) 
             sqqb_bbbtatann_EWp = sqqb_bbbtatann_EWp + t
              if (t .gt. 0d0 .and. .not. goodhel(ihel)) then
                  goodhel(ihel)=.true.
!                   write(*,*) ihel,t
              endif
          endif
      enddo
      sqqb_bbbtatann_EWp = sqqb_bbbtatann_EWp /  4d0
      end
c ============================ end function ====================================

c ============================== function ======================================     
c function generated by madgraph
c returns amplitude squared summed/avg over colors
c for the point in phase space p1,p2,p3,p4,p5,p6,p7,p8
c and helicity nhel(1),nhel(2)
c for process : q q -> b b l+ l- v_l v_l via (via A ,Z ,{Z'})
c ------------------------------- start  ---------------------------------------
      real*8 function qqb_bbbtatann_EWp(iq,gZpq,gZpt
     &                     ,Zpmass,Zpwidth,p1,p2,p3,p4,p5,p6,p7,p8,nhel)

c ----------------------------- declarations -----------------------------------
      implicit none

c constants
      integer    ngraphs ,nexternal
      parameter (ngraphs=7 ,nexternal=8)
      real*8     zero
      parameter (zero=0d0)

c arguments 
      integer iq ! incoming quark type (up/down)
      real*8 p1(0:3),p2(0:3),p3(0:3),p4(0:3),
     & p5(0:3),p6(0:3),p7(0:3),p8(0:3)  
      integer nhel(nexternal)
      real*8 gZpq(2,5),gZpt(2,5)
      real*8 Zpmass(5), Zpwidth(5)
          
c local variables 
      integer i,j
      integer jq
      parameter (jq=6) ! final state tops
      complex*16 amp_tmp
      complex*16 amp(ngraphs)
      complex*16 w1(6)  , w2(6)  , w3(6)  , w4(6)  , w5(6)
      complex*16 w6(6)  , w7(6)  , w8(6)  , w9(6)  , w10(6)
      complex*16 w11(6) , w12(6) , w13(6) , w14(6) , w15(6)
      real*8 gAq(2),gAt(2) ! coupling of a to q, t
      real*8 gZq(2),gZt(2) ! coupling of z to q, t
      real*8 gZpq_tmp(2),gZpt_tmp(2) ! necessary to pass 2d arrays
      
c global variables
      real*8         gW, gWWA, gWWZ
      common /coup1/ gW, gWWA, gWWZ
      real*8         gAl(2),gAu(2),gAd(2),gWf(2)
      common /coup2a/gAl,   gAu,   gAd,   gWf
      real*8         gZn(2),gZl(2),gZu(2),gZd(2),g1(2)
      common /coup2b/gZn,   gZl,   gZu,   gZd,   g1
      real*8         gWWh,gZZh,ghhh,gWWhh,gZZhh,ghhhh
      common /coup3/ gWWh,gZZh,ghhh,gWWhh,gZZhh,ghhhh
      complex*16     gh(2,12)
      common /coup4/ gh
      real*8         Wmass,Wwidth,Zmass,Zwidth
      common /vmass1/Wmass,Wwidth,Zmass,Zwidth
      real*8         Amass,Awidth,hmass,hwidth
      common /vmass2/Amass,Awidth,hmass,hwidth
      real*8            fmass(12), fwidth(12)
      common /fermions/ fmass,     fwidth
      integer     npoints
      common/stat/npoints
      integer       iQCD,iEW,iBSM
      common/igauge/iQCD,iEW,iBSM
      integer             iint
      common/interference/iint
c ---------------------------------- method ------------------------------------
c up/down type couplings (1-6: quarks, 7-12: leptons)

c initial down type
      if((iq.eq. 1).or.(iq.eq. 3).or.(iq.eq. 5).or.
     &   (iq.eq. 7).or.(iq.eq. 9).or.(iq.eq.11))then
        do i=1,2
          gaq(i)=gad(i)
          gzq(i)=gzd(i)
        end do
c initial up type
      else if((iq.eq. 2).or.(iq.eq. 4).or.(iq.eq. 6).or.
     &        (iq.eq. 8).or.(iq.eq.10).or.(iq.eq.12))then
        do i=1,2
          gaq(i)=gau(i)
          gzq(i)=gzu(i)
        end do
      end if
c final down type
      if((jq.eq. 1).or.(jq.eq. 3).or.(jq.eq. 5).or.
     &   (jq.eq. 7).or.(jq.eq. 9).or.(jq.eq.11))then
        do i=1,2
          gAt(i)=gad(i)
          gZt(i)=gzd(i)
        end do
c final down type
      else if((jq.eq. 2).or.(jq.eq. 4).or.(jq.eq. 6).or.
     &        (jq.eq. 8).or.(jq.eq.10).or.(jq.eq.12))then
        do i=1,2
          gAt(i)=gau(i)
          gZt(i)=gzu(i)
        end do
      end if

c initialise all the square matrix elements to zero
      do i=1,ngraphs,1
        amp(i)=0d0
      enddo


c wavefunctions
      call ixxxxx( p1  ,fmass(3  ) ,nhel(1  ) , 1 ,w1 )
      call oxxxxx( p2  ,fmass(3  ) ,nhel(2  ) ,-1 ,w2 )
      call oxxxxx( p3  ,fmass(12 ) ,nhel(3  ) , 1 ,w3 )
      call ixxxxx( p4  ,fmass(12 ) ,nhel(4  ) ,-1 ,w4 )
      call ixxxxx( p5  ,fmass(9  ) ,nhel(5  ) ,-1 ,w5 )
      call oxxxxx( p6  ,fmass(9  ) ,nhel(6  ) , 1 ,w6 )
      call oxxxxx( p7  ,fmass(10 ) ,nhel(7  ) , 1 ,w7 )
      call ixxxxx( p8  ,fmass(10 ) ,nhel(8  ) ,-1 ,w8 )

c W coupled to e- nu-bar vector current
      call jioxxx( w5 ,w7 ,gWf ,Wmass ,Wwidth ,w10 )
c W coupled to e+ nu vector current
      call jioxxx( w8 ,w6 ,gWf ,Wmass ,Wwidth ,w11 )
c t (off shell) coupled to b and W vector current 
      call fvoxxx( w3 ,w10 ,gWf ,fmass(11) ,fwidth(11) ,w12 )
c t-bar coupled to b-bar and W vector current 
      call fvixxx( w4 ,w11 ,gWf ,fmass(11) ,fwidth(11) ,w13 )

      if (iEW.eq.1)then
c A coupled to q q-bar vector current 
        call jioxxx( w1 ,w2 ,gAq ,Amass ,Awidth ,w9  )
c A diagram      
        call iovxxx( w13 ,w12 ,w9 ,gAt ,amp(1) )
c Z (off shell) coupled to q and q-bar vector current
        call jioxxx( w1 ,w2 ,gZq ,Zmass ,Zwidth ,w14 )
c Z diagram
        call iovxxx( w13 ,w12 ,w14 ,gZt ,amp(2) )
      else
        continue
      end if
c Z' diagrams
      if (iBSM.eq.1)then
        do i =1,5
          if (Zpmass(i).gt.0) then  
            do j=1,2
              gZpq_tmp(j)=gZpq(j,i)
              gZpt_tmp(j)=gZpq(j,i)
            end do     
            call jioxxx( w1 ,w2 ,gZpq_tmp ,Zpmass(i) ,Zpwidth(i) ,w15 )
            call iovxxx( w13 ,w12 ,w15   ,gZpt_tmp ,amp(2+i) )
          else
            continue
          end if 
        end do
      else
        continue  
      end if
      !write(*,*)'466'
c check amplitudes
!       do i=1,ngraphs,1
!         write(*,*)'M: ' ,i ,amp(i)
!       enddo

c total M*M for given helicity combination
      qqb_bbbtatann_EWp = 0.d0 
      amp_tmp = (0.d0,0.d0)
      if (iint.eq.0)then ! no interference
        do i=1,ngraphs
          qqb_bbbtatann_EWp = qqb_bbbtatann_EWp+amp(i)*conjg(amp(i))
        end do
      else if (iint.eq.1)then ! SM interference
        do i = 1, 2       
          amp_tmp = amp_tmp + amp(i)
        end do
        qqb_bbbtatann_EWp =qqb_bbbtatann_EWp+amp_tmp*conjg(amp_tmp)
        do i=3,ngraphs
          qqb_bbbtatann_EWp = qqb_bbbtatann_EWp+amp(i)*conjg(amp(i))
        end do
      else if (iint.eq.2)then ! full interference
        do i = 1, ngraphs       
          amp_tmp = amp_tmp + amp(i)
        end do
        qqb_bbbtatann_EWp =qqb_bbbtatann_EWp+amp_tmp*conjg(amp_tmp)
      else if (iint.eq.3)then ! interference only
        do i = 1, ngraphs       
          amp_tmp = amp_tmp + amp(i)
        end do
        qqb_bbbtatann_EWp =qqb_bbbtatann_EWp+amp_tmp*conjg(amp_tmp)
        do i=3,ngraphs
          qqb_bbbtatann_EWp = qqb_bbbtatann_EWp-amp(i)*conjg(amp(i))
        end do
      else
        write(*,*)'Error: interference flag not set.'
        stop
      end if
      
c print individual amplitudes
!       if (npoints.lt.10)then
!         do i=1,ngraphs
!           write(*,*)'M: ' ,i ,amp(i)
!         enddo
!       end if

      end
c ============================ end function ====================================
c ============================ end file ========================================