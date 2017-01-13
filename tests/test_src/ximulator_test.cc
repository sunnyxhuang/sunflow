//
// Created by Xin Sunny Huang on 10/22/16.
//

#include "gtest/gtest.h"
#include "src/events.h"

class XimulatorTest : public ::testing::Test {
 protected:
  virtual void TearDown() {
  }

  virtual void SetUp() {
    // TODO: replace ${YOUR_DIR_TO_SUNFLOW} with the absolute path to your
    // directory of the sunflow package,
    // e.g. "/Users/yourid/Documents/research/sunflow"
    // string TEST_DATA_DIR = "${YOUR_DIR_TO_SUNFLOW}/tests/test_data/";
    string TEST_DATA_DIR = "../test_data/";

    TRAFFIC_TRACE_FILE_NAME = TEST_DATA_DIR + "test_trace.txt";
    TRAFFIC_AUDIT_FILE_NAME = TEST_DATA_DIR  + "audit_traffic.txt";
    COMPTIME_AUDIT_FILE_NAME = TEST_DATA_DIR  +  "comp_time.txt";
    PORT_AUDIT_FILE_NAME = TEST_DATA_DIR  + "audit_port.txt";
    NET_TX_AUDIT_FILE_NAME = TEST_DATA_DIR  + "audit_net_tx.txt";
    CIRCUIT_AUDIT_FILE_NAME = TEST_DATA_DIR  + "audit_circuits.txt";

    ximulator_.reset(new Simulator());
  }
  std::unique_ptr<Simulator> ximulator_;
};

TEST_F(XimulatorTest, SunflowOnIntraCoflow) {
  ximulator_->InstallScheduler("sunflow");
  ximulator_->InstallTrafficGen("fb1by1");
  ximulator_->Run();
  EXPECT_NEAR(ximulator_->GetTotalCCT(), 3.141798, 1e-6);
}

TEST_F(XimulatorTest, SunflowOnInterCoflow) {
  ximulator_->InstallScheduler("sunflow");
  ximulator_->InstallTrafficGen("fb");
  ximulator_->Run();
  EXPECT_NEAR(ximulator_->GetTotalCCT(), 3.288448, 1e-6);
}


TEST_F(XimulatorTest, SolsticeOnIntraCoflow) {
  ximulator_->InstallScheduler("Solstice");
  ximulator_->InstallTrafficGen("fb1by1");
  ximulator_->Run();
  EXPECT_NEAR(ximulator_->GetTotalCCT(), 3.535875, 1e-6);
}

TEST_F(XimulatorTest, VarysOnInterCoflow) {
  ximulator_->InstallScheduler("varysImpl");
  ximulator_->InstallTrafficGen("fb");
  ximulator_->Run();
  EXPECT_NEAR(ximulator_->GetTotalCCT(), 2.710896, 1e-6);
}

TEST_F(XimulatorTest, AaloOnInterCoflow) {
  ximulator_->InstallScheduler("aaloImpl");
  ximulator_->InstallTrafficGen("fb");
  ximulator_->Run();
  EXPECT_NEAR(ximulator_->GetTotalCCT(), 2.928036, 1e-6);
}

TEST_F(XimulatorTest, BvNOnIntraCoflow) {
  ximulator_->InstallScheduler("BvN");
  ximulator_->InstallTrafficGen("fb1by1");
  ximulator_->Run();
  EXPECT_NEAR(ximulator_->GetTotalCCT(), 3.555875, 1e-6);
}

TEST_F(XimulatorTest, EdmondOnIntraCoflow) {
  ximulator_->InstallScheduler("Edmond");
  ximulator_->InstallTrafficGen("fb1by1");
  ximulator_->Run();
  EXPECT_NEAR(ximulator_->GetTotalCCT(), 4.481880, 1e-6);
}

TEST_F(XimulatorTest, TMSOnIntraCoflow) {
  ximulator_->InstallScheduler("TMS");
  ximulator_->InstallTrafficGen("fb1by1");
  ximulator_->Run();
  EXPECT_NEAR(ximulator_->GetTotalCCT(), 5.143361, 1e-6);
}

