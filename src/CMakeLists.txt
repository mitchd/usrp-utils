include_directories(${USRPutils_SOURCE_DIR}/include)

#Setup the programs
set(usrp_energy_SOURCES usrp-energy/usrp-energy.cpp)
set(usrp_recorder_SOURCES usrp-recorder/usrp-recorder.cpp)
set(usrp_sensor_SOURCES usrp-sensor/usrp-sensor.cpp common/fft_thread.cpp)
set(energycalculator_SOURCES energycalculator/energycalculator.cpp)
set(fftcompute_SOURCES fftcompute/fftcompute.cpp common/fft_thread.cpp)

add_executable(usrp_energy ${usrp_energy_SOURCES})
target_link_libraries(usrp_energy uhd m)

add_executable(energycalculator ${energycalculator_SOURCES})
target_link_libraries(energycalculator m)

add_executable(usrp_recorder ${usrp_recorder_SOURCES})
target_link_libraries(usrp_recorder uhd)

add_executable(usrp_sensor ${usrp_sensor_SOURCES})
target_link_libraries(usrp_sensor uhd m rt pthread fftw3f)

add_executable(fftcompute ${fftcompute_SOURCES})
target_link_libraries(fftcompute m rt pthread fftw3f)
