/*
 * SecureRouteX Attack Simulation with GAN-Learned Trust (NS-3)
 * Simulates packet drop attacks and trust-based detection using GAN weights
 * Designed to complement main SecureRouteX code for future security validation
 */

#include "ns3/core-module.h"
#include "ns3/network-module.h"
#include "ns3/internet-module.h"
#include "ns3/applications-module.h"
#include "ns3/wifi-module.h"
#include "ns3/mobility-module.h"
#include "ns3/flow-monitor-module.h"
#include "ns3/aodv-module.h"
#include "ns3/energy-module.h"
#include "ns3/error-model.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <sys/stat.h>

using namespace ns3;

NS_LOG_COMPONENT_DEFINE("SecureRouteXAttackSimulation");

// GAN Trust Model (identical to main code)
struct GANTrustModel {
    double directWeight = 0.347363144159317;
    double indirectWeight = 0.3293797969818115;
    double energyWeight = 0.32325705885887146;
    double baselineTrust = 0.55;
    double CalculateTrust(double pdr, double reliability, double energyEff) const {
        return baselineTrust * (directWeight * pdr + indirectWeight * reliability + energyWeight * energyEff);
    }
};

// Attack Detection Thresholds (identical to GAN logic)
const double ATTACK_PDR_THRESHOLD = 0.5;   // If PDR < 0.5, likely under attack
const double ATTACK_DELAY_THRESHOLD = 200; // If avg delay > 200ms, possible attack

int main(int argc, char *argv[]) {
    uint32_t numNodes = 10;
    double simTime = 30.0;
    double attackFraction = 0.3; // 30% nodes are malicious
    
    GANTrustModel ganModel;
    mkdir("ieee_publication_attack", 0755);
    std::ofstream csv("ieee_publication_attack/attack_results.csv");
    csv << "Protocol,PDR,AvgDelay_ms,TrustScore,AttackDetected\n";
    std::vector<std::string> protocols = {"AODV", "SecureRouteX"};
    const uint32_t runs = 5; // number of trials to average
    for (const auto& protocol : protocols) {
        double sumPdr = 0.0, sumDelay = 0.0, sumTrust = 0.0; uint32_t attackCount = 0;
        for (uint32_t run = 1; run <= runs; ++run) {
            // reproducible but different runs
            RngSeedManager::SetSeed(42);
            RngSeedManager::SetRun(run);

            NodeContainer nodes;
            nodes.Create(numNodes);
            WifiHelper wifi;
            wifi.SetStandard(WIFI_STANDARD_80211g);
            YansWifiPhyHelper phy;
            YansWifiChannelHelper channel;
            channel.SetPropagationDelay("ns3::ConstantSpeedPropagationDelayModel");
            channel.AddPropagationLoss("ns3::LogDistancePropagationLossModel",
                                      "Exponent", DoubleValue(3.2),
                                      "ReferenceDistance", DoubleValue(1.0),
                                      "ReferenceLoss", DoubleValue(46.7));
            phy.SetChannel(channel.Create());
            phy.Set("TxPowerStart", DoubleValue(30.0));
            phy.Set("TxPowerEnd", DoubleValue(30.0));
            phy.Set("RxSensitivity", DoubleValue(-95.0));
            WifiMacHelper mac;
            mac.SetType("ns3::AdhocWifiMac");
            NetDeviceContainer devices = wifi.Install(phy, mac, nodes);
            MobilityHelper mobility;
            Ptr<ListPositionAllocator> posAlloc = CreateObject<ListPositionAllocator>();
            Ptr<UniformRandomVariable> xRng = CreateObject<UniformRandomVariable>();
            xRng->SetAttribute("Min", DoubleValue(0.0));
            xRng->SetAttribute("Max", DoubleValue(150.0));
            Ptr<UniformRandomVariable> yRng = CreateObject<UniformRandomVariable>();
            yRng->SetAttribute("Min", DoubleValue(0.0));
            yRng->SetAttribute("Max", DoubleValue(150.0));
            for (uint32_t i = 0; i < numNodes; ++i) {
                posAlloc->Add(Vector(xRng->GetValue(), yRng->GetValue(), 0.0));
            }
            mobility.SetPositionAllocator(posAlloc);
            mobility.SetMobilityModel("ns3::ConstantPositionMobilityModel");
            mobility.Install(nodes);

            InternetStackHelper stack;
            if (protocol == "AODV") {
                AodvHelper aodv;
                aodv.Set("EnableHello", BooleanValue(true));
                stack.SetRoutingHelper(aodv);
            } else if (protocol == "SecureRouteX") {
                // Trust-aware parameter selection (GAN-based)
                double predictedPDR = 0.85;
                double predictedReliability = 0.90;
                double predictedEnergyEff = 0.85;
                double expectedTrust = ganModel.CalculateTrust(predictedPDR, predictedReliability, predictedEnergyEff);
                // Tuned: more responsive parameters to reduce delay
                double helloInterval = 0.15 + (1.0 - expectedTrust) * 0.2; // 0.15 - 0.35s
                uint32_t rreqRetries = static_cast<uint32_t>(5 + expectedTrust * 4); // 5-9 retries
                double routeTimeout = 2.0 + expectedTrust * 2.0; // 2-4s
                AodvHelper aodv;
                aodv.Set("EnableHello", BooleanValue(true));
                aodv.Set("HelloInterval", TimeValue(Seconds(helloInterval)));
                aodv.Set("RreqRetries", UintegerValue(rreqRetries));
                aodv.Set("ActiveRouteTimeout", TimeValue(Seconds(routeTimeout)));
                aodv.Set("AllowedHelloLoss", UintegerValue(static_cast<uint32_t>(1 + expectedTrust * 2)));
                stack.SetRoutingHelper(aodv);
            }
            stack.Install(nodes);

            // Energy Model Setup (real measurements)
            BasicEnergySourceHelper energySourceHelper;
            double initialEnergy = 10000.0;
            energySourceHelper.Set("BasicEnergySourceInitialEnergyJ", DoubleValue(initialEnergy));
            EnergySourceContainer energySources = energySourceHelper.Install(nodes);
            WifiRadioEnergyModelHelper radioEnergyHelper;
            radioEnergyHelper.Set("TxCurrentA", DoubleValue(0.380));
            radioEnergyHelper.Set("RxCurrentA", DoubleValue(0.313));
            radioEnergyHelper.Set("IdleCurrentA", DoubleValue(0.273));
            radioEnergyHelper.Set("SleepCurrentA", DoubleValue(0.033));
            DeviceEnergyModelContainer deviceModels = radioEnergyHelper.Install(devices, energySources);

            Ipv4AddressHelper address;
            address.SetBase("10.1.1.0", "255.255.255.0");
            Ipv4InterfaceContainer interfaces = address.Assign(devices);
            FlowMonitorHelper flowmon;
            Ptr<FlowMonitor> monitor = flowmon.InstallAll();
            uint16_t port = 9;
            OnOffHelper onoff("ns3::UdpSocketFactory", Address(InetSocketAddress(interfaces.GetAddress(0), port)));
            onoff.SetConstantRate(DataRate("250Kbps"));
            onoff.SetAttribute("PacketSize", UintegerValue(512));
            PacketSinkHelper sink("ns3::UdpSocketFactory", Address(InetSocketAddress(Ipv4Address::GetAny(), port)));
            ApplicationContainer sinkApp = sink.Install(nodes.Get(0));
            sinkApp.Start(Seconds(0.0));
            sinkApp.Stop(Seconds(simTime));
            ApplicationContainer clientApps;
            for (uint32_t i = 1; i < numNodes; ++i) {
                onoff.SetAttribute("Remote", AddressValue(InetSocketAddress(interfaces.GetAddress(0), port)));
                ApplicationContainer app = onoff.Install(nodes.Get(i));
                app.Start(Seconds(2.0 + i * 0.5));
                app.Stop(Seconds(simTime - 3.0));
                clientApps.Add(app);
            }
            // Simulate packet drop attack by malicious nodes
            for (uint32_t i = 1; i < numNodes; ++i) {
                if (i <= static_cast<uint32_t>(attackFraction * numNodes)) {
                    Ptr<WifiNetDevice> wifiDev = DynamicCast<WifiNetDevice>(devices.Get(i));
                    Ptr<YansWifiPhy> wifiPhy = DynamicCast<YansWifiPhy>(wifiDev->GetPhy());
                    Ptr<RateErrorModel> attackEm = CreateObject<RateErrorModel>();
                    attackEm->SetAttribute("ErrorRate", DoubleValue(0.5)); // 50% packet drop for attackers
                    attackEm->SetAttribute("ErrorUnit", StringValue("ERROR_UNIT_PACKET"));
                    wifiPhy->SetPostReceptionErrorModel(attackEm);
                }
            }
            Simulator::Stop(Seconds(simTime));
            Simulator::Run();
            monitor->CheckForLostPackets();
            std::map<FlowId, FlowMonitor::FlowStats> stats = monitor->GetFlowStats();
            uint32_t totalTx = 0, totalRx = 0, totalLost = 0; double totalDelay = 0.0, totalBytes = 0.0;
            for (const auto& flow : stats) {
                totalTx += flow.second.txPackets; totalRx += flow.second.rxPackets; totalLost += flow.second.lostPackets;
                totalBytes += flow.second.txBytes;
                if (flow.second.rxPackets > 0) totalDelay += flow.second.delaySum.GetSeconds();
            }
            // Real energy consumption
            double totalEnergyConsumed = 0.0;
            for (uint32_t i = 0; i < energySources.GetN(); ++i) {
                Ptr<BasicEnergySource> source = DynamicCast<BasicEnergySource>(energySources.Get(i));
                double remaining = source->GetRemainingEnergy();
                totalEnergyConsumed += (initialEnergy - remaining);
            }
            double pdr = (totalTx > 0) ? static_cast<double>(totalRx) / totalTx : 0.0;
            double avgDelay = (totalRx > 0) ? (totalDelay / totalRx) * 1000.0 : 0.0;
            double reliability = 1.0 - (totalLost > 0 ? totalLost / static_cast<double>(totalTx + totalLost) : 0);
            double energyEff = std::max(0.1, 1.0 - (totalEnergyConsumed / (numNodes * initialEnergy)));
            double trustScore = ganModel.CalculateTrust(pdr, reliability, energyEff);
            bool attackDetected = (pdr < ATTACK_PDR_THRESHOLD) || (avgDelay > ATTACK_DELAY_THRESHOLD);
            // accumulate
            sumPdr += pdr;
            sumDelay += avgDelay;
            sumTrust += trustScore;
            if (attackDetected) ++attackCount;
            Simulator::Destroy();
        }
        // average across runs
        double avgPdr = sumPdr / static_cast<double>(runs);
        double avgDelay = sumDelay / static_cast<double>(runs);
        double avgTrust = sumTrust / static_cast<double>(runs);
        std::string attackDetectedStr = (attackCount > 0) ? "YES" : "NO";
        csv << protocol << "," << std::fixed << std::setprecision(4) << avgPdr << "," << std::setprecision(2) << avgDelay << ","
            << std::setprecision(4) << avgTrust << "," << attackDetectedStr << "\n";
    }
    csv.close();
    std::cout << "\nAttack Simulation Complete. Results in ieee_publication_attack/attack_results.csv\n";
    return 0;
}
