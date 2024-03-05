$(OBJDIR)/frame.o           :frame.f FRAMELP;                        $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrlog_block.o   :mntrlog_block.f MNTRLOGD;               $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrlog.o         :mntrlog.f MNTRLOGD;                     $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrtmr_block.o   :mntrtmr_block.f MNTRLOGD MNTRTMRD;      $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/mntrtmr.o         :mntrtmr.f MNTRLOGD MNTRTMRD FRAMELP;    $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/rprm_block.o      :rprm_block.f RPRMD;                     $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/rprm.o            :rprm.f RPRMD FRAMELP;                   $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/io_tools_block.o  :io_tools_block.f IOTOOLD;               $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/io_tools.o        :io_tools.f IOTOOLD;                     $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/chkpoint.o        :chkpoint.f CHKPOINTD;                   $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/chkpt_mstp.o      :chkpt_mstp.f CHKPTMSTPD CHKPOINTD;      $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/map2D.o           :map2D.f MAP2D FRAMELP;                  $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/stat.o            :stat.f STATD MAP2D FRAMELP;             $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/stat_IO.o         :stat_IO.f STATD MAP2D FRAMELP;          $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/math_tools.o      :math_tools.f;                           $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/gSyEM.o           :gSyEM.f GSYEMD;                         $(FC) -c $(FL2) -I./ $< -o $@
$(OBJDIR)/io_trunc.o        :$S/3rd_party/io_trunc.f;                $(FC) -c $(FL2) -I./ $< -o $@
