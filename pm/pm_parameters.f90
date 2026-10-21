module pm_parameters

  logical,parameter :: part_memory=.true. ! Optimize particle memory distribution

  integer,parameter :: action_kick_only = 1
  integer,parameter :: action_kick_drift = 2

  integer,parameter :: PART_TYPE = 0
  integer,parameter :: STAR_TYPE = 1
  integer,parameter :: SINK_TYPE = 2
  integer,parameter :: TREE_TYPE = 3
  integer,parameter :: TRAC_TYPE = 4
  integer,parameter :: DUST_TYPE = 5
  integer,parameter :: CORE_TYPE = 6

end module pm_parameters
