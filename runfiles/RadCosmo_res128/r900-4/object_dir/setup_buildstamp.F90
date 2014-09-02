       subroutine setup_buildstamp (s_stamp_str, b_stamp_str, str_len)
       implicit none
       integer :: str_len
       character(len=str_len) :: s_stamp_str, b_stamp_str
       s_stamp_str = 'Tue Jun 3 13:30:59 2014'
       b_stamp_str = 'Tue Jun 3 13:31:10 2014'
       return
       end subroutine
      
       subroutine setup_systemInfo (system_str, str_len)
       integer :: str_len
       character(len=str_len) :: system_str
       system_str = 'Linux&
& r900-4.as.utexas.edu&
& 2.6.32-279.2.1.el6.x86_64&
& #1 SMP Thu Jul 5 21:08:58 EDT 2012&
& x86_64'
       return
       end subroutine
      
