import opacplot2.adapt
import opacplot2.opg_sesame
import os.path

def test_EosMergeGrids():
    # Check that the density grid satisfies lambda.
    ses_file = os.path.join(os.path.dirname(__file__), 'data/matr_009999.ses')
    op = opacplot2.opg_sesame.OpgSesame(ses_file, 
                                        opacplot2.opg_sesame.OpgSesame.SINGLE)
    eos_data = op.data[9999]
    merged_data = opacplot2.adapt.EosMergeGrids(eos_data,
                                          filter_temps=lambda x: x > 1,
                                          filter_dens=lambda x: x>1,
                                          intersect=['ele', 'ioncc'])
                                          # Setting intersect just in case
                                          # the defaults change & testsuite
                                          # doesn't.
    merged_data.run_testsuite(mode='full')